# © This script is modified from Eduard Fadeev @eddie1986
# NMDS - a non-parameteric monotonic relationship between the 
# dissimilarities in the samples matrix, and plots the location 
# of each site in a low-dimensional space (similar to principle component analysis).
#####################################
#setup
#####################################
setwd("~/Documents/...")
load("./data_import_v4v5.RData")
load("./data_import_v3v4.RData")

#####################################
#Load libraries
#####################################
library(ggplot2)
library(ggrepel)
library(phyloseq)
library(splitstackshape)
library(vegan)
library(ggdendro)
library(DESeq2)
################################################################################
# Centered Log-Ratio (CLR) Transformation
################################################################################
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

all_data.clr <- transform_sample_counts(all_data, fun = clr)

##################################### 
#Run "NMDS" using "bray"
#Run "PCoA" using "clr" or "euclidean"
#####################################
#Calculate distances

X1_nmds <- ordinate(physeq = all_data.clr, method = "NMDS", distance = "euclidean")
# X1_nmds <- ordinate(physeq = all_data, method = "NMDS", distance = "bray")
X1_nmds$stress # stress level of the plot if NMDS
#N.B. stress value below .2 is considered acceptable

# organize your parameters
sample_data(all_data)$Station<- as.factor(sample_data(all_data)$Station)
# sample_data(all_data.sub)$Substrate<- as.factor(sample_data(all_data.sub)$Substrate)
all_data@sam_data$Location <- factor(all_data@sam_data$Location, levels=c("Reference Site","Northern Seamount","Central Seamount","Karasik Seamount"))
all_data@sam_data$Depth <- factor(all_data@sam_data$Depth, levels=c("Surface","200 m","400 m", "600 m","Bathypelagic","Bottom Water"))


all_data@sam_data$Location <- factor(all_data@sam_data$Location, levels=c("Polaris Vent","Northern Seamount","Central Seamount","Seamount Saddle","Karasik Seamount","Southern Slope","S-Reference"))
all_data@sam_data$Layer <- factor(all_data@sam_data$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))
# all_data.sub_surf <- subset_samples(all_data.sub, Layer=="Surface (0-1cm)")

plot_ordination(
  physeq = all_data,
  ordination = X1_nmds,
  shape =  "Depth", #For water column
  # shape =  "Substrate", #For Sediment
  # shape =  "Layer", #For Sediment
  color = "Location"
  ) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Depth_m))+
  theme_classic() +
  # annotate("text", x = -50, y = -50, label = "Stress=0.03") + ##!!!!!!! IMPORTANT: CHECK STRESS LEVEL X1_nmds$stress
  annotate("text", x = -140, y = -100, label = "Stress=0.15") + ##!!!!!!! IMPORTANT: CHECK STRESS LEVEL X1_nmds$stress
  scale_x_continuous(limits = c(-150, 100))+
  scale_color_manual(values = c("Reference Site"="black","Northern Seamount"="#b4430f", "Central Seamount"="#E58601", "Karasik Seamount"="#0A5091"))+
  # scale_color_manual(values = c("Polaris Vent"="#ADADAD", "Northern Seamount"="#b4430f","Central Seamount"="#E58601","Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091","Southern Slope"="#1ebc97","S-Reference"="black"))+
  #shapescale sediment
  # scale_shape_manual(values = c(15,17,18)) +
  scale_shape_manual(values = c(16,17,18,25,8,15)) +
  # stat_ellipse(type = "norm", linetype = 2)+
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.text= element_text(size = 12),
          legend.position = "left",
          plot.title = element_text(size=12, lineheight=.8, face="bold", hjust = 0.5),
          panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
          plot.background = element_rect(color = "transparent")
          )
  
ggsave("./Figures/v3v4_NMDS.pdf", width = 10, height = 10,bg = "transparent")

# save.image(file="./NMDS_CommunityClustering.RData")

#See Statistics.R for Anosim and Permanova calculations