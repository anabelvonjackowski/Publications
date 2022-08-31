
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## NMDS
############################################################################################

# © This script is modified from Eduard Fadeev @eddie1986
# NMDS - a non-parameteric monotonic relationship between the dissimilarities in the samples matrix, and plots the location 
# of each site in a low-dimensional space (similar to principle component analysis).


#Load libraries

library(ggplot2)
library(ggrepel)
library(phyloseq)
library(splitstackshape)
library(vegan)
library(ggdendro)
library(DESeq2)

############################################################################################
#Run "NMDS" using "bray"
#Run "PCoA" using "clr" or "euclidean"
############################################################################################

# Calculate distances
# subset <- subset_samples(pseq.abs)

X1_nmds <- ordinate(
  physeq = pseq.abs, 
  method = "NMDS", 
  distance = "bray")

# stress level of the plot if NMDS
# NB - stress value below 0.2 is considered acceptable
X1_nmds$stress 

# organize your parameters
# sample_data(pseq.abs)$Station<- as.factor(sample_data(pseq.abs)$Station)
pseq.abs@sam_data$Depth <- factor(pseq.abs@sam_data$Depth, levels=c("Surface","DCM","BDCM", "100m"))
pseq.abs@sam_data$Station <- factor(pseq.abs@sam_data$Station, levels=c("N5","N4","N3", "S3","HG4","HG1/2","SV4","SV2"))

############################################################################################
## Plot
############################################################################################

plot_ordination(physeq = pseq.abs,
  ordination = X1_nmds,
  color =  "Depth") +
  geom_point(size = 6) +
  geom_text(aes(label = Station), nudge_y= -0.05, color = "black")+
  theme_classic() +
  annotate("text", x = -0.8, y = -1.0, label = "Stress=0.08") +
  # scale_color_manual(values = c("Summer" = "#31a354","Fall" = "black"))+
  scale_color_manual(values = c("#b34f96","#71ae5d","#6a65bd","#bc4949"))+
  scale_x_continuous(limits=c(-1.0,1.0))+
  scale_y_continuous(limits=c(-1.0,1.0))+
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.text= element_text(size = 12),
          legend.position = "left",
          plot.title = element_text(size=12, lineheight=.8, face="bold", hjust = 0.5),
          panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
          plot.background = element_rect(color = "transparent"))
  
ggsave("./Figures/NMDS_Depth_Location_Depth.pdf", width = 11, height = 10,bg = "transparent")
