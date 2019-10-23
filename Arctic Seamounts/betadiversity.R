# © This script is written by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
setwd("~/Documents/...")
load("./Data_Sed_v3v4_Rdata/NMDS_CommunityClustering.RData")
#####################################
#Load libraries
#####################################
library(phyloseq)
library(compositions)
library(ggplot2)
library(ggrepel)
library(grid)
library(igraph)
library(dplyr)
library(plyr)
library(reshape2)
library(gridExtra)
library(DESeq2)
library(gtable)
source("./multiplot.R")
source("./miseqR.R")
source("./plot_composition.R")
#####################################
#Subset the data set based on script 2_preprocessing
#####################################
# all_data.sub <- subset_samples(all_data.rare)
# all_data.sub<- subset_samples(all_data.rare.10K) #rarefied to 10K reads
# all_data.sub<- subset_samples(all_data.filt3) #filtrated by presence
# all_data.sub<- subset_samples(all_data.median) #normalized to median 
#####################################
#Plot top taxa
#####################################
all_data.sub.ra <- transform_sample_counts(all_data, function(x) 100 * x/sum(x))
all_data.sub <- subset_samples(all_data.sub.ra, Depth == "Surface")
all_data.sub <- subset_samples(all_data.sub.ra, Layer == "Surface (0-1cm)")

top_20_otus <- names(sort(taxa_sums(all_data.sub), TRUE)[1:20])
top_10_otus <- names(sort(taxa_sums(all_data.sub.ra), TRUE)[1:10])

X1_top.20 <- prune_taxa(top_20_otus, all_data.sub)
X1_top.10 <- prune_taxa(top_10_otus, all_data.sub.ra)

plot_bar(X1_top.20, y="Abundance", x="Location", fill="Class") 
plot_bar(X1_top.20, "Genus") 

##modified after https://github.com/joey711/phyloseq/issues/901

physeq2 = filter_taxa(all_data, function(x) mean(x) > 0.1, TRUE)
physeq2
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3

glom <- tax_glom(physeq3, taxrank = 'Class')
glom <- tax_glom(physeq3, taxrank = 'Genus', NArm=T)
tax_glom(physeq3, taxrank = "Order", NArm = FALSE)
glom # should list # taxa as # class

glom2 <- tax_glom(physeq3, taxrank = 'Genus', NArm=T)

data <- psmelt(glom) # create dataframe from phyloseq object
data2 = data

#simple way to rename phyla with < 1% abundance
data2$Class <- as.character(data2$Class) #convert to character
data2$Class[data$Abundance < 0.01] <- "Other Bacteria (<1%)"
data2$Class[data$Abundance < 0.01] <- "Other Archaea (<1%)"
data2$Class[data$Abundance < 0.02] <- "Other Archaea (<2%)"

#Count taxa you have 
Count = length(unique(data2$Class))
# Count = length(unique(data2$Genus))
Count

#Rename levels
data2$Location <- factor(data2$Location, levels=c("Polaris Vent","Northern Seamount", "Central Seamount","Seamount Saddle", "Karasik Seamount", "Southern Slope","Reference Site"))
# data2$Depth <- factor(data2$Depth, levels=c("Bottom Water","Bathypelagic", "600 m","400 m","200 m","Surface"))
data3$Depth <- factor(data3$Depth, levels=c("Surface", "200 m","400 m","600 m","Bathypelagic","Bottom Water"))
data2$Layer <- factor(data2$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))

data3=data2
#graph
ggplot(data=data3, aes(x = Location, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position="fill") +
  facet_grid(~Depth, switch = "y") + #y ist default
  ylab("Read Proportions") + 
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  xlab("")+
  guides(fill = guide_legend(ncol = 1))+ #1 or 2
  scale_fill_manual(values = my_color_class) + #See Color R Document
  # scale_x_discrete(position = "top")+ #if switch "x"
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 90,size = 12, hjust=0.95,vjust=0.2),
        # axis.text.x = element_text(angle = 90,size = 12, hjust=0,vjust=0.5), #if switch "x"
        axis.text.y = element_text(angle = 90,size = 12, vjust = 0,hjust = 0.5),
        strip.text = element_text(size = 12, angle = 90),
        legend.text = element_text(size = 8),
        legend.position= "right",
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1),
        strip.text.y = element_text(size = 12)
  )

ggsave("./figures/v4v5_Class_distribution_Location.png", width = 15, height = 10, dpi = 400)

save.image(file='./Data_Sed_v3v4_Rdata/betadiversity_barplots_v3v4.Rdata')
