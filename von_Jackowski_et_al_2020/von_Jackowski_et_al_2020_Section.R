#Plot the map using PlotSvalbard by @MikkoVihtakari
# https://github.com/MikkoVihtakari/PlotSvalbard

###################################
# Load Libraries
###################################
library(PlotSvalbard)
# install_github("MikkoVihtakari/PlotSvalbard", dependencies = TRUE)
# devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(extrafont)
# font_import(pattern="Palatino Linotype")

###################################
# Load Data
###################################
EW <-  data.merge %>%
  filter(Station %in% c("D1","D2","D3","D4_2","HG1","HG2","HG3", "HG4", "HG5","HG6","HG7", "HG8","HG9", "R2","SV1", "SV2", "SV3", "SV4"))

NS <- data.merge %>%
  filter(Station %in% c("N3", "N4", "N5", "R1", "R2", "R3", "S3","NSB_1"))

EWs.PS <- subset(EW, Cruise == "ARK32_PS114")
EWs.MSM <- subset(EW, Cruise == "MSM77")

###################################
# Plot
# jet.colors <- colorRampPalette(c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# plot.margin = Top, Right, Bottom, Left
###################################
section_plot(ARK.EW,
             z = "dHAA_umolL",
             x = "Longitude_degrees",
             y = "Depth_m", 
             sampling_indicator = "points",
             contour_label_cex = 2,
             interpolate = T)+
  xlab("Longitude (W/E)")+
  ylab("Depth (m)")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits = c(103,-5),trans = "reverse", breaks = c(0,25,50,100), expand = c(0,0))+ 
  scale_fill_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
                       # ,limits = c(0,5)
                       )+
  scale_color_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
                        # ,limits = c(0,5)
                        )
  #########
  # labs(color = bquote(atop('Chlorophyll-a','('*µg~L^-1*')')))+
  # labs(fill = bquote(atop('Chlorophyll-a','('*µg~L^-1*')')))+
  # labs(color = bquote(atop('Semi- \nlabile DOC','('*'%'~DOC*')')))+
  # labs(fill = bquote(atop('Semi- \nlabile DOC','('*'%'~DOC*')')))+
  # labs(color = bquote(atop('TEP Area','('*cm^2~L^-1*')')))+
  # labs(fill = bquote(atop('TEP Area','('*cm^2~L^-1*')')))+
  # labs(color = bquote(atop('Bacterial \nAbundance','('*10^5~cells~ml^-1*')')))+
  # labs(fill = bquote(atop('Bacterial \nAbundance','('*10^5~cells~ml^-1*')')))+
  # labs(color = bquote(atop('Bacterial \nProduction','('*µg~C~L^-1~d^-1*')')))+
  # labs(fill = bquote(atop('Bacterial \nProduction','('*µg~C~L^-1~d^-1*')')))+
  #########
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
        text=element_text(family="Palatino Linotype"),
        # axis.text.x=element_text(angle=0, hjust = 1),
        axis.title=element_text(size=40))

###################################
# Save
###################################
ggsave(filename="./Figures/PS114_Temperature.png", width = 13, height = 6, dpi= "retina")

Section <- ggarrange(PS.Chla, MSM.Chla,
                     align = "hv",
                     ncol=2, nrow=6)
ggsave("../Drafts/Latex/figures/Section_3.pdf",plot = Section, dpi= "retina", width = 35, height = 35)
