# Upload Data
setwd("~/Documents/PhD/My Papers/1_Seasonal 2018/R Work & Figures/")
#####################################
#Load Data
library(tidyr)
library(dplyr)
library(sp)
library(utils)
library(readxl)

# installation hack
library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("x/x")

#####################################
data_PS<- read.csv("./ARK_PS_2018.csv", sep = ";")

data_MSM<- read.csv("./ARK_MSM_2018.csv", sep = ";")

data.merge <- rbind(data_PS, data_MSM)

head(data.merge)
#####################################
#Subset
#####################################
EW <-  data %>%
    filter(Station %in% c("D1","D2","D3","D4_2","HG1","HG2","HG3", "HG4", "HG5","HG6","HG7", "HG8","HG9", "R2","SV1", "SV2", "SV3", "SV4"))
# or
EW <- subset(data.merge, Station != "N3" & Station != "N4" & Station != "N5" & Station != "S3" & 
                   Station != "R1" & Station != "R3" & Station != "R4" & Station != "R5")

NS <- data.merge %>%
    filter(Station %in% c("N3", "N4", "N5", "R1", "R2", "R3", "S3","NSB_1"))
# or
NS <- subset(data.merge, !(Station %in% EW$Station))


EWs.PS <- subset(EW, Cruise == "ARK32_PS114")

EWs.MSM <- subset(EW, Cruise == "MSM77")

NS.MSM <- subset(NS, Cruise == "MSM77")

#####################################
#further Structure manipulations
#####################################
data.merge$W_Layer <- factor(data.merge$W_Layer, levels=c("100m","bel.Chl","Chl.max","ab.Chl","Surface"))

data.merge$W_Layer <- factor(data.merge$W_Layer, levels=c("Surface","ab.Chl","Chl.max","bel.Chl","100m"))

data_subset<- na.omit(ARK.HG %>%
  dplyr::select(W_Layer,Station,Int_dCCHO_mmolm2) #%>%
    # subset(W_Layer == "Surface")
  )

#####################################
#count
#####################################
nrow(data_PS)

all.equal(dataset1,dataset2)
#####################################
#Write seperature CSV File
#####################################
write.csv(data_PS, file= "PS114.csv")