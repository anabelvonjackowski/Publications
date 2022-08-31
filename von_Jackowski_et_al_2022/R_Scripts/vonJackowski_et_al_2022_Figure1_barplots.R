############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## Stacked Barplot
############################################################################################

# Load packages
library(readxl)
library(tidyr) #for melt
library(ggplot2)

# Load Enrichment Dataframe

enrichment_processed <- as.data.frame(read_excel("./Enrichment.xlsx", sheet = "Enrichment"))
enrichment_processed$Parameter <-  factor(enrichment_processed$Parameter, levels=c("Bacteria (cells ml-1)","BP (umol C L-1 d-1)",
                                                                 "Synecoccus%POC","Synecoccus_cellsml","Picoplankton%POC","Picoplankton_cellsml","Phytoplankton (cells ml-1)",
                                                                 "PP (umol C L-1 d-1)","NDHAA%DOC","NDHAA_umol","EDHAA%DOC","EDHAA_umol",
                                                                 "AcDCCHO%DOC","AcDCCHO_umol","AmDCCHO%DOC","AmDCCHO_umol","NDCCHO%DOC","NDCCHO_umol",
                                                                 "DOC (umol L-1)","POC (umolL-1)"))

ggplot(enrichment_processed, aes(x=Enrichment, y= Parameter)) +
  geom_bar(stat = "identity",fill = "black",position = position_dodge(width=0.1)) +
  scale_x_continuous(trans = "log10", position = "bottom")+
  scale_y_discrete(position = "right") +
  geom_vline(xintercept=1)+
  scale_fill_manual(values = col.dccho)+
  guides(fill = guide_legend(nrow = 1))+
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
        # axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = .5),
        axis.title=element_text(size=12))

# Load Stacked Dataframe
Plankton_all <- na.omit(read.csv("./Github/Plankton_Stacked.csv"))

long_format <- gather(Plankton_all, key = "Plankton", value = "Relative",
                      c(15:18))
long_format$Plankton <- factor(long_format$Plankton, levels=c("Picoeukaryotes","Nanoeukaryotes","Cryptophytes","Synechococcus"))

long_format$Station <- factor(long_format$Station, levels = c("N5", "N4","N3","S3","HG4","HG1/2","SV4","SV2"))
long_format$Season <- factor(long_format$Season, levels=c("Summer","Fall"))
long_format$Layer <- factor(long_format$Layer, levels=c("Photic Zone","100 m"))
long_format$Depth <- factor(long_format$Depth_s, levels=c("Surface","DCM","Below DCM","100 m"))

ggplot(long_format, aes(x=Station, y= Relative, fill = Plankton)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Layer~Season, switch = "y")+
  # scale_x_continuous(trans = "reverse") +
  # scale_x_reverse()+
  scale_y_continuous(labels = scales::percent, expand = c(0,0),position = "right") +
  # scale_y_discrete(position = "right") +
  # coord_flip() +
  ylab("Abundance of Microbial Taxon")+
  scale_fill_manual(values = c("#b34f96","#71ae5d","#6a65bd","#bc4949"))+
  guides(fill = guide_legend(nrow = 1))+
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

ggsave("./Figures/Enrichment_Biogeochemistry.pdf", height = 9, width = 12)
