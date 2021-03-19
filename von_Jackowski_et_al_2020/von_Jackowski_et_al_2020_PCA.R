#PCA Variable Factor Map

#https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/

#####################################
#setup
#####################################
getwd()

#####################################
#Load libraries
#####################################
library(FactoMineR)
library(tidyverse)
library(devtools)
library(ggbiplot)
library(ggplot2)
library(plyr)
library(scales)
library(grid)

#####################################
#Load
#####################################
data_AA <- read.csv("./PCA/dHAA_PCA_Seasonal_Input.csv", header = TRUE,  sep = ",")
data_CHO <- read.csv("./PCA/dCCHO_PCA_Seasonal_Input.csv", header = TRUE,  sep = ",")

#####################################
#Run
#####################################
# exclude all data which are not AA [mol%], also exclude GABA
data <-data_AA[,c(1:13)]
ir.AA <- data_AA[,-c(14:20)] #additional data

data <-data_CHO[,c(1:10)]
ir.AA <- data_CHO[,-c(11:18)] #additional data

ir.AA$depth <- factor(ir.AA$depth, levels=c("Surface","above DCM","DCM","below DCM","100m"))
ir.AA$Timing <- factor(ir.AA$Timing, levels=c("Summer","Autumn"))


# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(data,center = TRUE, scale = TRUE) 
# print method
print(ir.pca)

ir.pca$rotation
y=as.data.frame(ir.pca$rotation) 

#export Results
# The print method returns the standard deviation of each of the four PCs,
# and their rotation (or loadings), which are the coefficients of the linear
# combinations of the continuous variables
write.csv(y,'output_PCA_Abb_Rotation.csv', row.names=T, col.names=T)


# more infos / output
# center and scale refers to respective mean and standard deviation of the variables that are used for normalization prior to implementing PCA
# The rotation measure provides the principal component loading. Each column of rotation matrix contains the principal component loading vector.
names(ir.pca)

ir.pca$x  #principal component score vector

x=as.data.frame(ir.pca$x) 

#export Results
write.csv(x, 'output_PCA_Abb_Scores.csv', row.names=T, col.names=T)


#####################################
#Plot
#####################################
#Simple Biplot
biplot(ir.pca, scale = 0)
dev.off()
# plot method
plot(ir.pca, type = "l")
# summary method
#The first row describe again the standard deviation associated with each PC.
#The second row shows the proportion of the variance in the data explained by each 
#component while the third row describe the cumulative proportion of explained variance. 
summary(ir.pca)


#Sophisticated Biplot
colnames(ir.AA)

ggbiplot(ir.pca, 
         choices = 1:2,
         varname.size = 5,
         var.axes = T,
         # groups = ir.AA$W_Layer, 
         # labels = ir.AA$W_Layer,
         obs.scale = 1, var.scale = 1
         ) +  
  # scale_color_manual(name="dCCHO [mol%]",values=c("navy","orange","#72FAA2","#A8FA72"))+
  geom_point(aes(shape=ir.AA$depth, fill=ir.AA$Timing), size=4)+
  scale_fill_manual(values=c("Summer"="cadetblue3","Autumn"="black"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  labs(shape="Depth", fill="Cruise")+
  scale_shape_manual(values=c(21,22,23,24,25))+
  ggtitle("dHAA [mol%]")+
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        # panel.background = element_rect(fill = "#faf0e6", color = NA), # bg of the plot
        panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        strip.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        # axis.text.x = element_text(angle = 90, hjust = 1, size=16),
        axis.title=element_text(size=12))

#####################################
#Save
#####################################
ggsave(file = "./PCA/dHAA_PCA.pdf", dpi = 300,height = 11, width = 15)
