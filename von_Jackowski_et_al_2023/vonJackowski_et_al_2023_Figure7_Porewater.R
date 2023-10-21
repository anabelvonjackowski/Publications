# © This script was written by Anabel von Jackowski @avonjackowski
# for the manuscript "Drivers of pelagic and benthic microbial communities on central Arctic seamounts"
#####################################
#setup
#####################################
rm(list = ls()) ##Clears entire workspace 
setwd("~/Documents/MarMic/Master_Thesis/Porewater_Profiles/")
scaleFUN <- function(x) sprintf("%.1f", x) # Sets decimal points for plots to 0.1
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
porewater <- read.csv("./PS101_station_Info_Porewater_2.csv", header = TRUE, check.names=FALSE)
porewater$Location <- factor(porewater$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle","Karasik Seamount", "Southern Slope", "S-Reference"))


#####################################
#Nitrate Profiles
#####################################
Nitrate <- ggplot(porewater, aes(y=Mean_Depth_cm, x= Nitrate))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(0,20),breaks=c(0,5,10,15,20)) + #x axis on top, set decimal, set scale
  ylab("Depth bsf (cm)")+
  xlab("µmol/L") +
  ggtitle("(c) Nitrate") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
Nitrate
#####################################
#Ammonium Profiles
#####################################
Ammonium <- ggplot(porewater, aes(y=Mean_Depth_cm, x= Ammonium))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(0, 8), breaks=c(0,2,4,6,8)) +
  ylab("Depth bsf (cm)")+
  xlab("µmol/L") +
  ggtitle("(e) Ammonium") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
Ammonium
#####################################
#Phosphate Profiles
#####################################
Phosphate <-ggplot(porewater, aes(y=Mean_Depth_cm, x= Phosphate))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(0, 1), breaks=c(0,.5,1.0)) +
  ylab("Depth bsf (cm)")+
  xlab("µmol/L") +
  ggtitle("(d) Phosphate") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
Phosphate
#####################################
#Silicate Profile
#####################################
Silicate <- ggplot(porewater, aes(y=Mean_Depth_cm, x= Silicate))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", breaks=seq(0, 160, by=40)) + #labels=scaleFUN, 
  ylab("Depth bsf (cm)")+
  xlab("µmol/L") +
  ggtitle("(f) Silicate") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
Silicate
#####################################
#Sulfate Profile
#####################################
Sulfate <- 
  ggplot(porewater, aes(y=Mean_Depth_cm, x= Sulfate))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(25, 35), breaks=c(25,30,35)) + #labels=scaleFUN, 
  ylab("Depth bsf (cm)")+
  xlab("mmol/L") +
  ggtitle("(g) Sulfate") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
Sulfate
#####################################
#Sulfide Profile
#####################################
Sulfide <- ggplot(porewater, aes(y=Mean_Depth_cm, x= Sulfide))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top", limits = c(0,3), breaks=c(0,1,2,3)) +
  ylab("Depth bsf (cm)")+
  xlab("µmol/L") +
  ggtitle("(h) Sulfide") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))

Sulfide
#####################################
#DIC
#####################################
DIC <- ggplot(porewater, aes(y=Mean_Depth_cm, x= c_DIC))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+ #reverse y axis
  scale_x_continuous(position = "top", limits = c(0, 3), breaks=c(0,1,2,3)) +
  ylab("Depth bsf (cm)")+ #y axis label
  xlab("mmol/L") + #x axis label
  ggtitle("(a) DIC") + #title
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
DIC
#####################################
#TA Profile
#####################################
TA <-ggplot(porewater, aes(y=Mean_Depth_cm, x= c_alcalinity))+
  facet_wrap(~Location, strip.position="bottom", nrow = 1) + #, margins = TRUE) + 
  theme_bw()+
  geom_path(size=1,aes(color =  Station))+ #defines path color to be by location
  geom_point(aes(color = Station), stat="identity", shape = 4, size = 2)+
  scale_color_manual(values = c("PS101/140"="#C27D38","PS101/187"="grey80",
                                "PS101/211"="#0B775E", "PS101/101"="#0A5091","PS101/123"="#5ca3e6",
                                "PS101/151"="#ffcc66", "PS101/218"="#1ebc97", "PS101/063"="black")) +
  scale_y_reverse()+
  scale_x_continuous(position = "top",limits = c(0, 3), breaks=c(0,1,2,3)) +
  ylab("Depth bsf (cm)")+
  xlab("mmol/L") +
  ggtitle("(b) Total Alkalinity") +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=0),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))
TA
#####################################
#Save - Compile files into one large image
#####################################
# ggsave("./Sulfide_Porewater_20180216.pdf", plot = Sulfide)

porewater.all <- grid.arrange(DIC,TA,Nitrate,Phosphate,Ammonium,Silicate,Sulfate,Sulfide,
                              ncol=2, nrow=4)
ggsave("./Plots/Supplementary_Porewater.pdf", dpi = 300, plot = porewater.all,  units = "mm",
       width = 297, height = 210)
