# © This script was written by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
rm(list = ls()) ##Clears entire workspace 
setwd("~/Documents/MarMic/Master_Thesis/Porewater_Profiles/")
#####################################
#Load libraries
#####################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(ggpubr)
library(devtools)
#####################################
#Load Data
#####################################

sediment <- read.csv("./PS101_station_Info_Sediment.csv", header = TRUE, check.names=FALSE)

#####################################
#Plot
#####################################

ggplot(sediment)+
  theme_bw()+
  facet_grid(~Location, switch="both") + # strip.position="bottom", margins = TRUE) +
  geom_point(aes(y=Depth_cm, x= Phaeopigment,fill = Location), stat="identity", shape = 4, size = 3)+
  # geom_path(aes(y=Depth_cm, x= Phaeopigment, color =  Location), size = 1)+
  scale_color_manual(values = c("Rift Valley"="#ADADAD", "Reference Site"="black",
                                "Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091",
                                "Central Mount"="#E58601","Northern Mount"="#b4430f", "Southern Slope"="#1ebc97"))+
  scale_fill_manual(values = c("Rift Valley"="#ADADAD", "Reference Site"="black",
                               "Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091",
                               "Central Mount"="#E58601","Northern Mount"="#b4430f", "Southern Slope"="#1ebc97"))+
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0,1.5), breaks=c(0,.5,1,1.5)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/ml")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("Phaeopigment")

ggsave("./Plots/Chlorophyll_Profile.pdf", plot = profile,width = 7, height = 12)
CPE.all <- ggarrange(CPE_Plot, Phae_plot, Chla_plot, Cont_Plot, TOC_Plot, TON_Plot, C_N_plot,ncol=2, nrow=4)
ggsave("./Plots/Supplementary Figure 4.png",plot = CPE.all, width = 25, height = 15)
