#Map based on PlotSvalbard

######################
#Load Packages
######################
library(readxl)
library("PlotSvalbard")



library(measurements)
library(ggplot2)
library(ggrepel)
library(oce)

######################

######################
#Set up data
######################
data <- data.merge

# convert degree-min to decimal degree
# data$latDD <- gsub('°', ' ', data$lat)
# data$lonDD <- gsub('°', ' ', data$long)
data$lonDD <- gsub('W', '-', data$Longitude_degrees)
data$lonDD <- gsub('E', '', data$Longitude_degrees)
data$latDD <- data$Latitude_degrees

data$lonDD <- as.numeric(data$lonDD)

sta <- data.frame(lon = data$lonDD, lat = data$latDD, st = data$Station)
sta <- unique(sta)

head(sta)
str(sta)

######################
#Plot the map using PlotSvalbard by @MikkoVihtakari
# install_github("MikkoVihtakari/PlotSvalbard", dependencies = TRUE)
#https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html
######################
stas <- transform_coord(sta, lon = "lon", lat = "lat", bind = T) # for svarlbard projection

# barentssea mapping is "+init=epsg:32633"
#limits = c(1*10^6, -1*10^6, 0, -2*10^6) for a stenographic map or  limits = c(-30, 30, 75, 85) for Greenland Sea
#"arctic60" for a stenographic map or "svalbard" for Greenland Sea
# bathymetry = T or F

basesta <- basemap("barentssea", limits = c(-10, 20, 75, 82), 
                   currents = T,
                   current.size = "scaled", 
                   legends = c(T, F),
                   bathymetry = F,
                   land.col = "grey90",
                   land.border.col = "black")
basesta

######################
#Plot
######################
 basesta+
  geom_point(data = stas.t, aes(x=lon.utm, y=lat.utm), pch=19, size=1, color = "black")


######################
#Plot
######################
ggsave(filename="./Figures/Map.pdf")