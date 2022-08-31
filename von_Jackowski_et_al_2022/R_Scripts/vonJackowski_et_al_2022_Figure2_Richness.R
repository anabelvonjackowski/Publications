
##############################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## rarefaction and diversity analyses
##############################################################################

# set working directory; save and load
# setwd("/AWI_MPI/FRAM/collaborations/PEBCAO/Rstats")
# load("Pebcao.Rdata")

library(iNEXT) # install.packages("iNEXT")
library(olsrr)
library(cowplot)
library(phyloseq)
library(ggplot2)

############################################################################################
### Calculate Results
############################################################################################
iNEXT <- otu_table(ASV, taxa_are_rows = F)

iNEXT <- iNEXT(as.data.frame(otu_table(iNEXT)), q=c(0),
               datatype="abundance", conf = 0.95, nboot = 100)

rich <- iNEXT$AsyEst[iNEXT$AsyEst$Diversity == "Species richness",]
shan <- iNEXT$AsyEst[iNEXT$AsyEst$Diversity == "Shannon diversity",]
simp <- iNEXT$AsyEst[iNEXT$AsyEst$Diversity == "Simpson diversity",]

AlphaDiv <- data.frame(Fram_ID = iNEXT$DataInfo$site)
AlphaDiv$Depth <- ENV$Depth 
AlphaDiv$Depth <- factor(AlphaDiv$Depth, levels = c("Surface", "Photic Zone","Base Photic Zone","100m"))
AlphaDiv$Layer <- ENV$Layer
AlphaDiv$watermass <- ENV$watermass 
AlphaDiv$Season = ENV$Season
AlphaDiv$Station = ENV$Station
AlphaDiv$Station <- factor(AlphaDiv$Station, levels = c("N5", "N4","N3","S3","HG4","HG1/2","SV4","SV2"))
AlphaDiv$Sample_sum = iNEXT$DataInfo$n
AlphaDiv$Observed = iNEXT$DataInfo$S.obs
AlphaDiv$Richness = rich$Observed
AlphaDiv$Richness.cov = rich$Observed/rich$Estimator
AlphaDiv$Shannon = shan$Observed
AlphaDiv$Shannon.est = shan$Estimator
AlphaDiv$Simpson = simp$Observed
AlphaDiv$Simpson.est = simp$Estimator
AlphaDiv$Sam.comp = 100*iNEXT$DataInfo$SC

############################################################################################  
### Export and Plot
############################################################################################
write.table(AlphaDiv, file="AlphaDiversity.txt", sep="\t", row.names=F)


# AlphaDiv.photic <- subset(AlphaDiv, Layer == "Photic")
# AlphaDiv.100m <- subset(AlphaDiv, Layer == "100m")

# Diversity 
ggplot(data = AlphaDiv, aes(x = Station, y = Richness)) +
  geom_boxplot() +
  # scale_y_continuous(limits = c(500,2000))+
  geom_point(aes(color = Depth), size = 6) +
  # geom_point(data = AlphaDiv.100m,aes(x = Station, y = Shannon,color = Depth), size = 6) +
  scale_color_manual(values = c("Surface"="#b34f96","DCM"="#71ae5d","Below DCM"="#6a65bd","100m"="#bc4949"))+
  # scale_color_manual(values=col.season) +
  facet_grid(.~Season)+ #, switch = "y"
  # scale_y_continuous(position = "right")+
  # geom_text(aes(label = Station))+
  # ylab("Richness")+
  theme_bw() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = .5),
        axis.title=element_text(size=12))

############################################################################################
### Save
############################################################################################
# Single Graphic 
ggsave("./Figures/Alpha_Richness_depth.pdf",width = 8, height = 5,  plot = Shannon, bg = "transparent")

# Graphic Array
library(ggpubr)
SummaryPlot <- ggarrange(Richness, Shannon, Simpson,align = "v",ncol=1, nrow=3, common.legend = T, legend="bottom")
ggsave("./Figures/Alpha_all.pdf",plot = SummaryPlot, width = 6, height = 9)

############################################################################################
### Statistical Analysis
############################################################################################
anova <- aov(Richness  ~ Season, data = AlphaDiv)
summary(anova)
TukeyHSD(anova)

anova <- aov(Richness ~ Layer, data = subset(AlphaDiv,Season=="Summer"))
summary(anova)
TukeyHSD(anova)

## PERMANOVA
anova <- aov(Simpson ~ Season*Depth, data = AlphaDiv)
summary(anova)
TukeyHSD(anova)

## ANOSIM: continuous variables
library(vegan)
adonis(Shannon ~ Season*Depth, method = "bray", data = AlphaDiv)

adonis(Simpson ~ Season*Depth, method = "bray",
       # data = subset(AlphaDiv, Season=="Fall")
       data = AlphaDiv)


############################################################################################
### Delete temporary datasets
rm(rich, simp, shan, cover)
rm(list=ls(pattern = "cover.point.*|cover.line.*|rarefac.point.*|rarefac.line.*"))
############################################################################################

save.image("Pebcao.Rdata")
