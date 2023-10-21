# © This script was written by Anabel von Jackowski @avonjackowski
# for the manuscript "Drivers of pelagic and benthic microbial communities on central Arctic seamounts"

#####################################
#setup
#####################################
# rm(list = ls()) ##Clears entire workspace 
# setwd("~/Documents/...)

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
sediment <- read.csv("./PS101_Sediment.csv", header = TRUE, check.names=FALSE)
sediment$Contribution <- sediment$Contribution*100
sediment$Location <- factor(sediment$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle","Karasik Seamount", "Southern Slope", "S-Reference"))
sediment <- na.omit(sediment)
Chlorophyll <-subset(sediment, select = c(Location, Replicate, Mean_Depth_cm, Chlorophyll))
Chlorophyll$Location <- factor(Chlorophyll$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle", "Karasik Seamount", "Southern Slope", "S-Reference"))

#####################################
#Profile CPE
#####################################
CPE.profile <-
  ggplot(sediment)+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= CPE,color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(y=Mean_Depth_cm, x= CPE, color =  Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 2), breaks=c(0,1,2)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/ml")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("CPE")

#####################################
#Profile Phaeopigment
#####################################
Phae.profile <-
  ggplot(sediment)+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= Phaeopigment,color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(y=Mean_Depth_cm, x= Phaeopigment, color =  Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 2), breaks=c(0,1,2)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/ml")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("CPE")

#####################################
#Profile Chlorophyll
#####################################
Chla.profile <-
  ggplot(sediment)+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= Chlorophyll,color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(y=Mean_Depth_cm, x= Chlorophyll, color =  Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 0.5), breaks=c(0,0.25,0.5)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/ml")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("Chlorophyll")

#####################################
# Chlorophyll Contribution
#####################################
Contribution <- ggplot(sediment, aes(y=Mean_Depth_cm, x= Contribution))+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= Contribution, color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(color = Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top") +
  ylab("Depth bsf (cm)") +
  xlab("%")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  ggtitle("Chlorophyll Contribution (%)")

#####################################
#Load Data
#####################################
sediment <- read.csv("./PS101_Sediment.csv", header = TRUE, check.names=FALSE)
sediment$Location <- factor(sediment$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle", "Karasik Seamount", "Southern Slope", "S-Reference"))
TOC_TON <- subset(sediment, select = -c(Chlorophyll,Phaeopigment, CPE, Contribution))
TOC_TON <- na.omit(TOC_TON)

#####################################
# TOC Concentration
#####################################
TOC <- ggplot(TOC_TON)+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= TOC,color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(y=Mean_Depth_cm, x= TOC, color = Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 13), breaks=c(0,4,8,12)) +
  # scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 2), breaks=c(0,1,2)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/mg")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("Total Organic Carbon (TOC)")

#####################################
# TON Concentration
#####################################
TON <- ggplot(TOC_TON)+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= TON,color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(y=Mean_Depth_cm, x= TON, color = Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", labels=scaleFUN, limits = c(0, 2), breaks=c(0,1,2)) +
  ylab("Depth bsf (cm)") +
  xlab("µg/mg")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12))+
  ggtitle("Total Organic Nitrogen (TON)")

#####################################
#Ratio C:N
#####################################
CtoN <- ggplot(TOC_TON, aes(y=Mean_Depth_cm, x= C_N))+
  theme_bw()+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  geom_point(aes(y=Mean_Depth_cm, x= C_N, color = Replicate), stat="identity", shape = 4, size = 4)+
  geom_path(aes(color =  Replicate), size = 1)+
  scale_color_manual(values = c("black","#5ca3e6","#ffcc66","#899DA4","#C27D38","forestgreen")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(0, 16), breaks=c(0,4,8,12,16)) +
  ylab("Depth bsf (cm)") +
  xlab("C:N")+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 12)) +
  ggtitle("C:N ratio")

#####################################
#Save Files
#####################################
# ggsave("./Plots/Chlorophyll_Profile.pdf", plot = profile,width = 7, height = 12)

CPE.all <- grid.arrange(CPE.profile,Phae.profile,Chla.profile,Contribution,nrow=2)
ggsave("./Plots/Supplementary_Chlorophyll.pdf",plot = CPE.all,
       dpi = 300, units = "mm",width = 297, height = 105)

TOC_TON_plot <- grid.arrange(TOC, TON, CtoN, nrow=2)
ggsave("./Plots/Supplementary_OM.pdf", dpi = 300, plot = TOC_TON_plot,  units = "mm",
       width = 297, height = 105)