
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## Section Plot
############################################################################################

library(PlotSvalbard)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(readxl)
library(measurements)

###################################
###  Plot
###################################

data_load <- read.csv(file = "./metadata.csv", header = TRUE,  sep = ",")

data_EW <-  data_load %>%filter(Station %in% c("HG1","HG4", "SV2", "SV4"))

data_NS <- data_load %>% filter(Station %in% c("N3", "N4", "N5", "S3"))

EW.PS <- subset(EW, Cruise == "ARK32_PS114")
EW.MSM <- subset(EW, Cruise == "MSM77")

NS.PS <- subset(NS, Cruise == "ARK32_PS114")
NS.MSM <- subset(NS, Cruise == "MSM77")


section_plot(NS.MSM,
             z = "Temperature_C",
             x = "Longitude_degrees",
             x = "Latitude_degrees",
             y = "Depth_m", 
             sampling_indicator = "points",contour_label_cex = 5,interpolate = T)+
  ylab("Depth (m)")+
  xlab("Longitude (W/E)")+
  scale_x_continuous(expand = c(0,0), limits=c(4,10))+
  xlab("Latitude (N)")+
  scale_x_continuous(expand = c(0,0), limits = c(78.5,80))+
  scale_y_continuous(limits = c(103,-5),trans = "reverse", breaks = c(0,25,50,100), expand = c(0,0))+ 
  scale_fill_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),limits = c(0,7.1))+
  scale_color_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),limits = c(0,7.1))+
  labs(color = bquote(atop('Temperature','('*C*')')))+
  labs(fill = bquote(atop('Temperature','('*C*')')))+
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
      plot.margin = unit(c(0,8,0,3), "cm"),
      # plot.title = element_text(hjust = -0.1,vjust = 0,size = 20),
      panel.background = element_rect(fill = "grey50", color = "white", size = 2),
      # legend.position = c(1.25, 0.51),
      legend.title = element_text(size=40), 
      # legend.title = element_text(size=25, hjust = 6), 
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      legend.background = element_rect(fill = "transparent", color = "transparent"), #legend bg
      legend.text = element_text(size=40),
      legend.key.size = unit(1.7, "cm"),
      axis.text=element_text(size=40),
      # axis.text.x=element_text(angle=0, hjust = 1),
      axis.title=element_text(size=40))

###################################
# Save
###################################
ggsave(filename="x.pdf", width = 13, height = 6, dpi= "retina")

Section <- ggarrange(x,x,x,x,x,x,x,
                     align = "hv",ncol=2, nrow=3)

ggsave("xx.pdf",plot = Section, dpi= "retina")
