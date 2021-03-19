# © This script was written by Anabel von Jackowski @avonjackowski
# See: https://rpubs.com/JanpuHou/278558

#####################################
#setup
#####################################
rm(list = ls()) ##Clears entire workspace 

scaleFUN <- function(x) sprintf("%.1f", x) # Sets decimal points for plots to 0.1

#####################################
#Load libraries
#####################################
library(ggplot2)
library(FactoMineR)
library(factoextra)

#####################################
#Load Data
#####################################
data_PS <- read.csv("ARK_PS_2018.csv", header = TRUE,  sep = ",")
data_MSM <- read.csv("ARK_MSM_2018.csv", header = TRUE,  sep = ",")

#####################################
#Cluster Analysis
#####################################

data.cluster <- subset(data.PS, select=c(Temperature_C,Salinity_PSU,Phytoplankton_cells_ml,Bacteria_cells_ml,HNA_ml,LNA_ml,
                                        Bac_Prod1_nmolLeuL.1h.1,Bac_Prod1_microgCl.1d.1,
                                        Cell_Growth_Rate_d.1,Turnover_d.1))
rownames(data.cluster) <- data.s$Stationbook
rownames(data.cluster) <- make.names(data.s[,1], unique = TRUE) #station
rownames(data.cluster) <- make.names(data.s[,4], unique = TRUE) #depth
rownames(data.cluster) <- make.names(data.s[,5], unique = TRUE) #Current

df <- scale(data.cluster)
km.res <- kmeans(na.omit((data.cluster)),4, nstart = 25)
aggregate(df, by=list(cluster=km.res$cluster), mean)

fviz_cluster(km.res, data = df,
             # palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot"
)
ggsave("./PS114_Cluster.pdf")

#####################################
#cluster dendrogram
# see: https://rpkgs.datanovia.com/factoextra/reference/fviz_dend.html
#####################################
# Hierarchical clustering
#Example
res <- hcut(USArrests, k = 4, stand = TRUE)
#Mine
res <- hcut(data.cluster, k = 4, stand = TRUE)

# Visualize
fviz_dend(res, rect = TRUE, cex = 1,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))
ggsave("./PS114_Dendrogram.pdf")
