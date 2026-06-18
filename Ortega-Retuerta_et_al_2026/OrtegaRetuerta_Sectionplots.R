################################################################################
# Project: OdiSea Section Visualization
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Updated: June 2026
# Description: Plot the map/transect profiles using PlotSvalbard by @MikkoVihtakari
# Repository: https://github.com/MikkoVihtakari/PlotSvalbard
################################################################################

###################################
# Load Libraries
###################################
library(readxl)       # For reading environmental profiles out of Excel sheets
library(PlotSvalbard) # For high-level oceanographic ocean section plots and spatial mapping
# install_github("MikkoVihtakari/PlotSvalbard", dependencies = TRUE)
# devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
library(ggplot2)      # For modular graphics layouts and adjustments
library(dplyr)        # For data manipulation and pipelined tracking workflows
library(ggpubr)       # For arranging multiple grid matrices or panel configurations
library(extrafont)    # For custom typefaces and external font installation support
# font_import(pattern="Palatino Linotype")
library(scales)       # For internal scale mapping breaks and palette visualizers

# Color Palette Resource Reference: https://bids.github.io/colormap/
show_col(viridis_pal()(7))
show_col(viridis_pal(option = "plasma")(7))

# Oceanographic Color Management Resource References:
# https://matplotlib.org/cmocean/
# https://cran.r-project.org/web/packages/cmocean/vignettes/cmocean.html
library(cmocean)
# Palette Options Validation Check Error Notes: 
# ‘name’ must be one of “algae”, “amp”, “balance”, “curl”, “deep”, “delta”, “dense”, 
# “diff”, “gray”, “haline”, “ice”, “matter”, “oxy”, “phase”, “rain”, “solar”, “speed”,
# “tarn”, “tempo”, “thermal”, “topo”, “turbid”
image(volcano, col = cmocean('algae')(7)) 
list(volcano, col = cmocean('algae')(7)) 
# "#D7F9D0FF" "#A1D494FF" "#64B463FF" "#129450FF" "#126E45FF" "#1A482FFF" "#122414FF"

###################################
# Load & Clean Ocean Data
###################################
# Isolate profile data matching the target cruise year "2021"
Profile <- subset(as.data.frame(read_excel('./mola_data_profile.xlsx', col_names = T)),
                  Year == "2021")

# Format columns into categorical factor levels and strictly numerical elements
Profile$Year <- as.factor(Profile$Year)
Profile$Salinity_psu <- as.numeric(Profile$Salinity_psu)
Profile$Depth <- as.numeric(Profile$Depth)
str(Profile)

###################################
# Plot Transect Profile
# Setup and Range Diagnostics:
# jet.colors <- colorRampPalette(c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# plot.margin = Top, Right, Bottom, Left
# colnames(Meta_subset)
###################################
# Diagnose sample vector boundaries before plotting
min(na.omit(Profile$Bacterial_105cellsmL))
max(na.omit(Profile$Bacterial_105cellsmL))

###################################
# Generate Interpolated Oceanographic Section Plot
section_plot(Profile,
             z = "Bacterial_105cellsmL",    # Response variable for contour colors
             x = "JulianDay",               # X-axis time element
             y = "Depth",                   # Y-axis depth element
             sampling_indicator = "points", # Overlay actual physical bottle samples
             contour = 5,                   # Interval spacing for overlay isolines
             contour_color = "white",       # Isoline color choice
             contour_label_cex = 1.2,       # Size threshold for numeric contour labels
             interpolate = T)+              # Enable grid matrix spatial interpolation
  xlab("Julian Day")+
  ylab("Depth (m)")+
  
  # Configure time horizons matching specific monthly transition days
  scale_x_continuous(expand = c(0,0), limits = c(60,244), breaks = c(60, 91, 121, 152, 182, 213, 244))+
  
  # Reverse Y-axis scales to place zero (ocean surface) at the top border boundary
  scale_y_continuous(expand = c(0,0), limits = c(500,0),trans = "reverse", breaks = c(0,100, 200, 300, 400 , 500))+
  
  # [Commented Out] Alternative standard, viridis, chl, or oxy color profiles
  # scale_fill_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))+
  # scale_fill_gradientn(colours = c("#042333FF", "#2C3396FF","#744992FF", "#B05F82FF", "#EB7958FF", "#FBB33DFF" ,"#E8FA5BFF"),
  # scale_fill_gradientn(colours = c("#440154FF","#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF"),
  # scale_fill_gradientn(colours = c("white", "#D7F9D0FF", "#64B463FF", "#129450FF", "#126E45FF", "#1A482FFF","#122414FF"),#chl
  # scale_fill_gradientn(colours = c("#400505FF", "#850A0BFF", "#6F6F6EFF", "#9A9A99FF", "#CBCAC9FF", "#EBF44CFF","#DDAF19FF"),#oxy
  
  # Fill Palette Settings: Apply custom 'deep' cmocean palette structure
  scale_fill_gradientn(colours = c("#3F396DFF","#3E6495FF","#488E9EFF", "#5DBAA4FF", "#A4DEA6FF", "white", "#FDFECCFF", "#FDE725FF"),#deep
                       limits = c(0,11))+
  
  # [Commented Out] Alternative point/line gradient profiles
  # scale_color_gradientn(colours = c("#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))+
  # scale_color_gradientn(colours = c("#042333FF", "#2C3396FF", "#744992FF", "#B05F82FF", "#EB7958FF", "#FBB33DFF", "#E8FA5BFF"),
  # scale_color_gradientn(colours = c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF"),
  # scale_color_gradientn(colours = c("white", "#D7F9D0FF", "#64B463FF", "#129450FF", "#126E45FF", "#1A482FFF", "#122414FF"), #chl
  # scale_color_gradientn(colours = c("#400505FF", "#850A0BFF", "#6F6F6EFF", "#9A9A99FF", "#CBCAC9FF", "#EBF44CFF","#DDAF19FF"), #oxy
  
  # Color Boundary Settings: Matching structural mapping color scale limits 
  scale_color_gradientn(colours = c("#3F396DFF","#3E6495FF","#488E9EFF", "#5DBAA4FF", "#A4DEA6FF", "white", "#FDFECCFF", "#FDE725FF"),#deep
                        limits = c(0,11))+
  # labs(color = bquote(atop('Chlorophyll','('*ugL*')')))+
  # labs(fill = bquote(atop('Chlorophyll','('*ugL*')')))+
  
  # Theme layout properties configured for seamless, transparent background rendering
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"), # bg of the plot
        # plot.margin = unit(c(0,8,0,3), "cm"),
        # plot.title = element_text(hjust = -0.1,vjust = 0,size = 20),
        panel.background = element_rect(fill = "transparent", color = "transparent", size = 2),
        # legend.position = c(1.25, 0.51),
        legend.title = element_text(size=12), 
        # legend.title = element_text(size=25, hjust = 6), 
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent", color = "transparent"), #legend bg
        legend.text = element_text(size=12),
        legend.key.size = unit(1, "cm"),
        axis.text=element_text(size=12),
        # text=element_text(family="Palatino Linotype"),
        # axis.text.x=element_text(angle=0, hjust = 1),
        axis.title=element_text(size=12))

###################################
# Save Figures
###################################
# ggsave(filename="./Figures_Section/MOLA_Section.pdf", width = 87.5, height = 70, units = "mm", dpi= "retina") #3.2
