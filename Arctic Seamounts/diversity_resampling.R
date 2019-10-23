# © This script was written by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
setwd("~/Documents/...")
load("./Data_Sed_v3v4_Rdata/NMDS_CommunityClustering.RData")
#####################################
#Load libraries
#####################################
source("~/Documents/MarMic/Master_Thesis/16S/scripts/SubSampleNGS_mod.R")
#####################################
#Calculate Alpha WITH replication
#http://www.evolution.unibas.ch/walser/bacteria_community_analysis/2015-02-10_MBM_tutorial_combined.pdf
#####################################
data = read.csv("./OTU_table.csv", header = TRUE, row.names = 1)
OTU_number = rowSums(t(data)) #total number of sequences for your samples that you are using
min(OTU_number) #smallest dataset
# calculating Alpha diversity indices based on n random subsampling runs to e.g. the minimum library size (sub)
Alpha <- SubSampleNGS(data, 100, min(OTU_number))
#####################################
#Run Diversities
#####################################
#Richness - count of species
Data01=data
Data01[Data01>0]<-1
nOTU = apply(Data01,2,sum)
#with resampling
nOTU_r <- Alpha[["nOTU"]]
nOTU_r_mean = as.dist(apply(nOTU_r,2,mean))
nOTU_r_sd = apply(nOTU_r,2,sd)

# Richness - estimates diversity from abundance data (importance of rare OTUs) 
# - assumes that the number of obervations for a taxa has a Poisson distribution
# and    corrects    for    variance;    
chao1 <- Alpha[["chao1"]]
chao1_mean = apply(chao1,2,mean) ###try to accumulate it into one table
chao1_sd = apply(chao1,2,sd)
#save.image("UnivV4V5_Chao1.Rdata")

inverseSimpson = invS=diversity(t(data), "inv")
#resampling
invS_r <- Alpha[["invS"]]
invS_r_mean = apply(invS_r,2,mean)
invS_r_sd = apply(invS_r,2,sd)
#save.image("UnivV4V5_invSimp.Rdata")
#shannon abundance and evenness
shannon <- diversity(t(data), "shannon")
#resampling
shan_r <- Alpha[["shannon"]]
shan_r_mean = apply(shan_r,2,mean)
shan_r_sd = apply(shan_r,2,sd)
#save.image("UnivV4V5_Shannon.Rdata")

##Do exponential shannon
shannon.exp <- exp(diversity(t(data), "shannon"))
#####################################
#mantel test
#https://sites.google.com/site/mb3gustame/hypothesis-tests/the-mantel-test
#build disimilarity matrix using vegdist
#h0=number of OTU
#h1=exponential shannon 
#h2=inverse simpon
#####################################
mantel.nOTU <- mantel(vegdist(nOTU, method = "euclidian"), vegdist(nOTU_r_mean, method = "euclidian"))
#V4V5 = Mantel statistic r:  0.98
#V4V5 = Significance: 0.001 
#water V3V4 = Mantel statistic r: -0.008964 
#water V3V4 = Significance: 0.575
#sed V3V4 = Mantel statistic r: 0.09583 
#sed V3V4 = Significance: 0.003
mantel.invS <- mantel(vegdist(inverseSimpson, method = "euclidian"), vegdist(invS_r_mean, method = "euclidian"))
#V4V5 = Mantel statistic r:     1 
#V4V5 = Significance: 0.001 
#water V3V4 = Mantel statistic r:     1  
#water V3V4 = Significance: 0.001 
#sed V3V4 = Mantel statistic r:     0.9914
#sed V3V4 = Significance: 0.001 
mantel.shannon <- mantel(vegdist(shannon, method = "euclidian"), vegdist(shan_r_mean, method = "euclidian"))
#V4V5 = Mantel statistic r:     1 
#V4V5 = Significance: 0.001 
#water V3V4 = Mantel statistic r: 0.9996  
#water V3V4 = Significance: 0.001
#sed V3V4 = Mantel statistic r: 0.9864  
#sed V3V4 = Significance: 0.001
#####################################
##Summaries
#####################################
all_div = cbind(nOTU, nOTU_r_mean, nOTU_r_sd, chao1_mean, chao1_sd, inverseSimpson, invS_r_mean, invS_r_sd, shannon, shan_r_mean, shan_r_sd)
data2 = cbind(nOTU, chao1_mean, chao1_sd, inverseSimpson, shannon, shannon.exp)
write.csv(all_div, file = "./Watercolumn_Diversityind_all_v3v4.csv")
all = list(nOTU, chao1_mean, chao1_sd, inverseSimpson, shannon)
