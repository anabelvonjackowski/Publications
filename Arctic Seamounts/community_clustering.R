# © This script is modified from Eduard Fadeev @eddie1986
# NMDS - a non-parameteric monotonic relationship between the 
# dissimilarities in the samples matrix, and plots the location 
# of each site in a low-dimensional space (similar to principle component analysis).
#####################################
#setup
#####################################
setwd("~/Documents/...")
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
##################################### 
#Run "NMDS" using "bray"
#Run "PCoA" using "clr" or "euclidean"
#####################################
#Calculate distances
all_data.sub <- subset_samples(all_data)

X1_nmds <- ordinate(
  physeq = all_data.sub, 
  method = "NMDS", 
  distance = "bray"
)
X1_nmds$stress # stress level of the plot if NMDS
#N.B. stress value below .2 is considered acceptable

# organize your parameters
sample_data(all_data.sub)$Station<- as.factor(sample_data(all_data.sub)$Station)
sample_data(all_data.sub)$Substrate<- as.factor(sample_data(all_data.sub)$Substrate)
all_data.sub@sam_data$Location <- factor(all_data.sub@sam_data$Location, levels=c("Polaris Vent","Northern Seamount",
                                                                          "Central Seamount","Seamount Saddle",
                                                                          "Karasik Seamount","Southern Slope","Reference Site"))

all_data.sub@sam_data$Depth <- factor(all_data.sub@sam_data$Depth, levels=c("Surface","200 m","400 m", "600 m","Bathypelagic","Bottom Water"))
all_data.sub@sam_data$Layer <- factor(all_data.sub@sam_data$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))
all_data.sub_surf <- subset_samples(all_data.sub, Layer=="Surface (0-1cm)")

plot_ordination(
  physeq = all_data.sub,
  ordination = X1_nmds,
  # shape =  "Depth", #For water column
  # shape =  "Substrate", #For Sediment
  shape =  "Layer", #For Sediment
  color = "Substrate"
  ) +
  geom_point(size = 4) +
  #geom_text_repel(aes(label = Station))+
  theme_classic() +
  # annotate("text", x = 1.5, y = -1.5, label = "Stress=0.02") + ##!!!!!!! IMPORTANT: CHECK STRESS LEVEL X1_nmds$stress
  annotate("text", x = 1.0, y = -1.5, label = "Stress=0.1") + ##!!!!!!! IMPORTANT: CHECK STRESS LEVEL X1_nmds$stress
  # scale_color_manual(values = c("Polaris Vent"="#ADADAD", "Reference Site"="black",
  #                                "Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091","Karasik Seamount N"="#0A5091",
  #                                "Central Seamount"="#E58601","Northern Seamount"="#b4430f", "Southern Slope"="#1ebc97"))+
  #shapescale sediment
  # scale_shape_manual(values = c(19,18,15,17)) +
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
  
ggsave("./figures/NMDS_Depth_Location.pdf", width = 10, height = 10,bg = "transparent")

save.image(file="./NMDS_CommunityClustering.RData")

#See Statistics.R for Anosim and Permanova calculations