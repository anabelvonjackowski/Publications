# © This script was written by Anabel von Jackowski @avonjackowski
#####################################
#setup
# setwd("~/Documents/MarMic/Master_Thesis/Porewater_Profiles/")
library(vegan)
#####################################

#####################################
#Import Data
#####################################
sediment <- read.csv("./PS101_station_Info_Sediment.csv", header = TRUE) 
data<- subset(sediment, select = c(Nitrate,Ammonium,Phosphate,Silicate,Sulfate,TA, DIC,
                                   Location_g,Location, Substrate_g,Substrate, Depth_cm,Depth_g))
data$Depth_g <- factor(data$Depth_g, levels=c("Surface", "Subsurface", "Deep Layer"))
data <- na.omit(data)

#####################################
# Simple PCA
#####################################
data("iris")
head(iris)
#plot the pairwise scatterplot
pairs(data,lower.panel = NULL)
#Perform PCA
profile.pca <- princomp(data, cor = TRUE, scores = TRUE)
#cor: a logical value. If TRUE, the data will be centered and scaled before the analysis
#scores: a logical value. If TRUE, the coordinates on each principal component are calculated
summary(profile.pca)
plot(profile.pca)
#Run biplot
biplot(profile.pca) 
#expand
biplot(profile.pca, expand=10, xlim=c(-0.30, 0.0), ylim=c(-0.1, 0.1))
#Run Loadings. These are coefficients in linear combination (with the eigenvectors?) predicting a variable by the (standardized) components
profile.pca$loadings
#Get Scores. These is the transformed dataset
profile.pca$scores

#####################################
#more complex PCA
library("FactoMineR")
library("ggplot2")
library("factoextra")
#####################################
data.pca <- PCA(data[,1:7], scale.unit = T,ncp = 5, graph = FALSE)

PCA_Location <- fviz_pca_biplot(data.pca,
                axes = c(1, 2), 
                geom.ind = "point", #leaving this out causes default geom = c("point", "text"),
                geom.var = c("arrow", "text"),
                col.ind = "white",
                fill.ind = data$Location,
                # alpha.ind=.5,
                # habillage = data.merged$W_Type, #coloring the observations by groups
                mean.point = FALSE,
                pointshape = 21,
                pointsize = 4,
                col.var = "black", #Arrow Color
                # fill.var = "white",
                gradient.cols = NULL,
                # gradient.cols = c("blue", "yellow", "red"),
                label = "all", 
                invisible = "none",
                repel = FALSE, 
                palette = NULL, #palette = my_color # palette = "jco"
                # addEllipses = FALSE,
                # ellipse.type = "confidence",
                title = ""
) +
  labs(fill = "Location", color = "Contribution")+
  ggtitle("")+
  scale_fill_manual(values = c("Rift Valley"="#ADADAD", "Abyssal "="black","Seamount"="#0A5091","Slope"="#1ebc97",
                              "Reference Site"="black","Seamount Saddle"="#0B775E", "Karasik Seamount"="#0A5091",
                                              "Central Seamount"="#E58601","Northern Seamount"="#b4430f", "Southern Slope"="#1ebc97",
                              "Sediment"="black", "Sediment and Basaltic Pebbles"="#ADADAD", "Sponge + Tubeworm Debris"="#E58601", "Tubeworm Debris"="#b4430f",
                              "Surface"="#E58601","Subsurface"="#ADADAD","Deep Layer"="black"))+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        text = element_text(size=12),
        strip.text = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title=element_text(size=12))

#####################################
#save
#####################################
PCA.all <- ggarrange(PCA_Location_g,PCA_Location,PCA_Substrate,PCA_Layer,ncol=2, nrow=2, common.legend = F, legend="bottom")
# ggsave("./PCA_Sediments_Porewater_Location_g.pdf")
ggsave("./Plots/Supplementary Figure 5.png",plot = PCA.all,width = 25, height = 15,dpi = 300)