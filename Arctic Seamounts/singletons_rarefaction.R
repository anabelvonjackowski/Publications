# © This script is written by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
setwd("~/Documents/...")
rm(list = ls()) ##Clears entire workspace 
#####################################
#Load libraries
#####################################
library(phyloseq)
library(ggplot2)
library(plyr)
library(reshape2) #for melt function
library(ggpubr) #for 2 plots
source("../scripts/rarefaction_curve.R")
#####################################
#Double check whether there are singletons in the dataset
#####################################
data = read.csv("./OTU_table.csv", header = TRUE, row.names = 1)
data_no_singletons = data[rowSums(data)!=1,] #double check whether there are singletons. If this and command abouve are equal, there are none
identical(data, data_no_singletons)
#####################################
#count
#####################################
OTU_number = rowSums(t(data)) #total number of sequences for your samples
min(OTU_number) #smallest dataset
number <- as.data.frame(OTU_number)
number$Sample <- rownames(number)
#####################################
## contribution of the singletons
# merge file - OTU number --> make percentage (sum of total database)
# weight of the singletons on the total database 
#####################################
nSeqs_merged = read.csv("nSeqs_V4V5_all.csv", header = 1)
#This require some manual cut and pastes as this script was written at the very beginning of the thesis
#In the CSV, just add a column that referres to the sample code "X.."

singletons <- merge(nSeqs_merged, number, by = "Sample")
singletons$Total <- (singletons$Merged + singletons$OTU_number)
singletons$Singletons <- (singletons$Merged - singletons$OTU_number)
singletons$Contribution <- ((singletons$Singletons/singletons$Total)*100)
Singleton_Contribution <- singletons

save.image(file="./Data_Watercolumn_2_Rdata/singletons.RData") #water column
save.image(file="./Data_Sed_v3v4_Rdata/singletons.RData") #sediment
#####################################
#Plot rarefaction curves (adjusted from Bela @and3k )
#####################################
set.seed(42)

sample_sums(all_data)

#Calculate rarefaction curve
#possible indices are: 'Observed', 'Chao1', 'ACE', 'Shannon', 'Simpson', 'InvSimpson, 'Fisher'

rarefaction_curve_data <- calculate_rarefaction_curves(all_data, c('Observed', 'Chao1', 'Shannon','InvSimpson'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))

summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), 
                                        summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(all_data)), 
                                                by.x = 'Sample', by.y = 'row.names')
# Rename facets to correct names
# indices_names <- c(`Observed` = "Richness",`Chao1` = "Chao 1",`se.chao1` = "Chao 1 Standard Error",`Shannon` = "Shannon-Wiener",`InvSimpson` = "Inverse Simpson")

#specify order on plot
rarefaction_curve_data_summary_verbose$Location <- factor(rarefaction_curve_data_summary_verbose$Location, levels=c("Polaris Vent","Northern Seamount", "Central Seamount", "Seamount Saddle","Karasik Seamount",    "Southern Slope","Reference Site"))
#Depth + Layer
rarefaction_curve_data_summary_verbose$Depth.y <- factor(rarefaction_curve_data_summary_verbose$Depth.y, levels=c("Surface", "200 m","400 m","600 m", "Bottom Water", "Bathypelagic"))
all_data@sam_data$Layer <- factor(all_data@sam_data$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))
#Convert "Station" from integer to factor
rarefaction_curve_data_summary_verbose$Station <- as.factor(rarefaction_curve_data_summary_verbose$Station)


rarefraction_plot <- ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    # x = Depth,
    x = Depth.x,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Location,
    # shape = Depth.y, #for water column
    shape = Layer, #for sediment
    group = Sample)) + 
  geom_line(size =1) + 
  geom_pointrange(size = .3) + 
  facet_wrap(facets = ~ Measure, scales = 'free') +
  ylab("Alpha Diversity Mean") +
  xlab("Sequences") +
  scale_color_manual(values = c("Polaris Vent"="#ADADAD", "Reference Site"="black",
                                "Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091",
                                "Central Seamount"="#E58601","Northern Seamount"="#b4430f", "Southern Slope"="#1ebc97")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        strip.text = element_text(size = 12),
        legend.text= element_text(size = 12),
        plot.title = element_text(size=12, lineheight=.8, face="bold", hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))

ggsave("./figures/rarefaction_v3v4.pdf", width = 12, height = 10)

alpha.divs <- ggarrange(richness.plot,rarefraction_plot ,ncol=1, nrow=2)
ggsave("./figures/Alphadiversity_V3v4.pdf",plot = alpha.divs, width = 8, height = 10)

