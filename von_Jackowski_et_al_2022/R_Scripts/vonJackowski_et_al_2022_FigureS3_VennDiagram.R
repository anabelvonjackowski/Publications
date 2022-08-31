
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## Venn Diagram
############################################################################################

library(venn)
library(ggpolypath)
library(ggplot2)
library(tibble)
library(gtools)
library(dplyr)

OTU_dataframe <- data.frame(ASV[,mixedsort(names(ASV))])
colnames(OTU_dataframe) = ENV$Fram_ID_A
OTU_dataframe_t <- t(OTU_dataframe)

OTU_dataframe <- as.data.frame(OTU_dataframe_t, add_rownames=T) %>%
  as_tibble(rownames = "ID")
# View(OTU_dataframe)

##############################################################################
###  Subset Dataframes
##############################################################################
Summer_All <- filter(OTU_dataframe, ID=="PS_HG4_1" | ID=="PS_HG4_2"| ID=="PS_HG4_3" | ID=="PS_HG4_4"| 
                                  ID=="PS_HG1_2_1" |ID=="PS_HG1_2_2" |  ID=="PS_HG1_2_3" |ID=="PS_HG1_2_4" |
                                  ID == "PS_SV4_1"|  ID == "PS_SV4_2"|ID == "PS_SV4_3"|  ID == "PS_SV4_4"|
                                  ID == "PS_SV2_1"|ID == "PS_SV2_2"| ID == "PS_SV2_3"|ID == "PS_SV2_4"|
                                  ID == "PS_N5_1"|ID == "PS_N5_2"|  ID == "PS_N5_3"|ID == "PS_N5_4"| 
                                  ID == "PS_N4_1" |ID == "PS_N4_2" | ID == "PS_N4_3" |ID == "PS_N4_4" |
                                  ID == "PS_N3_1" |ID == "PS_N3_2" |ID == "PS_N3_3" |ID == "PS_N3_4" |
                                  ID == "PS_S3_1" | ID == "PS_S3_2" |ID == "PS_S3_3" | ID == "PS_S3_4")
Summer_All <- t(Summer_All[,2:length(Summer_All)])
Summer_Photic <- filter(OTU_dataframe, ID=="PS_HG4_1" | ID=="PS_HG4_2"| ID=="PS_HG4_3"| 
                           ID=="PS_HG1_2_1"| ID=="PS_HG1_2_2" | ID=="PS_HG1_2_3" |  
                           ID == "PS_SV4_1"|  ID == "PS_SV4_2"| ID == "PS_SV4_3"|
                           ID == "PS_SV2_1"| ID == "PS_SV2_2"| ID == "PS_SV2_3"|
                           ID == "PS_N5_1"| ID == "PS_N5_2"| ID == "PS_N5_3"| 
                           ID == "PS_N4_1" | ID == "PS_N4_2" | ID == "PS_N4_3" |
                           ID == "PS_N3_1" | ID == "PS_N3_2" | ID == "PS_N3_3" |
                           ID == "PS_S3_1" | ID == "PS_S3_2" | ID == "PS_S3_3" )
# row.names(Summer_Surface_Photic) <- Summer_Surface_Photic$ID
Summer_Photic <- t(Summer_Photic[,2:length(Summer_Photic)])

Summer_100m <- filter(OTU_dataframe, ID=="PS_HG4_4" |
                        ID=="PS_HG1_2_4" |
                        ID=="PS_S3_4" | 
                        ID == "PS_SV4_4" |
                        ID=="PS_SV2_4" | 
                        ID=="PS_N5_4" |
                        ID=="PS_N4_4" |
                        ID=="PS_N3_4")
Summer_100m <- t(Summer_100m[,2:length(Summer_100m)])


Fall_All <- filter(OTU_dataframe, ID=="MSM_HG4_1" | ID=="MSM_HG4_2"|  ID=="MSM_HG4_3" |
                                ID=="MSM_HG1_1" |ID=="MSM_HG1_2" |  ID=="MSM_HG1_3" |ID=="MSM_HG1_4" | 
                                ID == "MSM_SV4_1"|  ID == "MSM_SV4_2"|ID == "MSM_SV4_3"|  ID == "MSM_SV4_4"|
                                ID == "MSM_SV2_1"|ID == "MSM_SV2_2"|ID == "MSM_SV2_3"|ID == "MSM_SV2_4"|
                                ID == "MSM_N5_1"|ID == "MSM_N5_2"| ID == "MSM_N5_3"|ID == "MSM_N5_4"| 
                                ID == "MSM_N4_1" |ID == "MSM_N4_2" |ID == "MSM_N4_3" |ID == "MSM_N4_4" |
                                ID == "MSM_N3_1" |ID == "MSM_N3_2" |ID == "MSM_N3_3" |ID == "MSM_N3_4" |
                                ID == "MSM_S3_1" | ID == "MSM_S3_2"|  ID == "MSM_S3_3" | ID == "MSM_S3_4")
Fall_All <- t(Fall_All[,2:length(Fall_All)])
Fall_Photic <- filter(OTU_dataframe, ID=="MSM_HG4_1" | ID=="MSM_HG4_2"| 
                         ID=="MSM_HG1_1" | ID=="MSM_HG1_2" | ID=="MSM_HG1_3" | 
                         ID == "MSM_SV4_1"| ID == "MSM_SV4_2"| ID == "MSM_SV4_3"|
                         ID == "MSM_SV2_1"| ID == "MSM_SV2_2"| ID == "MSM_SV2_3"|
                         ID == "MSM_N5_1"| ID == "MSM_N5_2"| ID == "MSM_N5_3"| 
                         ID == "MSM_N4_1" | ID == "MSM_N4_2" | ID == "MSM_N4_3" |
                         ID == "MSM_N3_1" | ID == "MSM_N3_2" |ID == "MSM_N3_3" |
                         ID == "MSM_S3_1" | ID == "MSM_S3_2" | ID == "MSM_S3_3")
Fall_Photic <- t(Fall_Photic[,2:length(Fall_Photic)])

Fall_100m <- filter(OTU_dataframe, ID=="MSM_HG4_3"| 
                      ID=="MSM_HG1_4" | 
                      ID == "MSM_SV4_4"|
                      ID == "MSM_SV2_4"| 
                      ID == "MSM_N5_4"| 
                      ID == "MSM_N4_4" | 
                      ID == "MSM_N3_4" | 
                      ID == "MSM_S3_4")
Fall_100m <- t(Fall_100m[,2:length(Fall_100m)])

# Individuals
Summer_Surface <- filter(OTU_dataframe, ID=="PS_HG4_1"|ID=="PS_HG1_2_1" |
                           ID == "PS_SV4_1"|ID == "PS_SV2_1"|ID == "PS_N5_1"|
                           ID == "PS_N4_1" |ID == "PS_N3_1" |ID == "PS_S3_1")
Summer_Surface <- t(Summer_Surface[,2:length(Summer_Surface)])
Summer_Photic <- filter(OTU_dataframe, ID=="PS_HG4_2"|ID=="PS_HG1_2_2" | 
                          ID == "PS_SV4_2"|ID == "PS_SV2_2"|ID == "PS_N5_2"|
                          ID == "PS_N4_2" |ID == "PS_N3_2" |ID == "PS_S3_2")
Summer_Photic <- t(Summer_Photic[,2:length(Summer_Photic)])
Summer_BasePhotic <- filter(OTU_dataframe, ID=="PS_HG4_3"|ID=="PS_HG1_2_3" | 
                              ID == "PS_SV4_3"|ID == "PS_SV2_3"|ID == "PS_N5_3"|
                              ID == "PS_N4_3" |ID == "PS_N3_3" |ID == "PS_S3_3")
Summer_BasePhotic <- t(Summer_BasePhotic[,2:length(Summer_BasePhotic)])
Summer_100m <- filter(OTU_dataframe, ID=="PS_HG4_4"|ID=="PS_HG1_2_4" | 
                              ID == "PS_SV4_4"|ID == "PS_SV2_4"|ID == "PS_N5_4"|
                              ID == "PS_N4_4" |ID == "PS_N3_4" |ID == "PS_S3_4")
Summer_100m <- t(Summer_100m[,2:length(Summer_100m)])

Fall_Surface <- filter(OTU_dataframe, ID=="MSM_HG4_1"|ID=="MSM_HG1_2_1" |
                         ID == "MSM_SV4_1"|ID == "MSM_SV2_1"|ID == "MSM_N5_1"|
                         ID == "MSM_N4_1" |ID == "MSM_N3_1" |ID == "MSM_S3_1")
Fall_Surface <- t(Fall_Surface[,2:length(Fall_Surface)])
Fall_Photic <- filter(OTU_dataframe, ID=="MSM_HG4_2"|ID=="MSM_HG1_2_2" | 
                        ID == "MSM_SV4_2"|ID == "MSM_SV2_2"|ID == "MSM_N5_2"|
                        ID == "MSM_N4_2" |ID == "MSM_N3_2" |ID == "MSM_S3_2")
Fall_Photic <- t(Fall_Photic[,2:length(Fall_Photic)])
Fall_BasePhotic <- filter(OTU_dataframe, ID=="MSM_HG4_3"|ID=="MSM_HG1_2_3" | 
                            ID == "MSM_SV4_3"|ID == "MSM_SV2_3"|ID == "MSM_N5_3"|
                            ID == "MSM_N4_3" |ID == "MSM_N3_3" |ID == "MSM_S3_3")
Fall_BasePhotic <- t(Fall_BasePhotic[,2:length(Fall_BasePhotic)])
Fall_100m <- filter(OTU_dataframe, ID=="MSM_HG4_4"|ID=="MSM_HG1_2_4" | 
                      ID == "MSM_SV4_4"|ID == "MSM_SV2_4"|ID == "MSM_N5_4"|
                      ID == "MSM_N4_4" |ID == "MSM_N3_4" |ID == "MSM_S3_4")
Fall_100m <- t(Fall_100m[,2:length(Fall_100m)])

##############################################################################
###  Venn Diagram
##############################################################################
# Summer vs Fall
comb <- data.frame(cbind(rowSums(Summer_All), rowSums(Fall_All))>0)
names(comb)<-c("Summer","Fall")

# Summer Photic and Fall Photic vs Summer 100 and Fall 100
comb <- data.frame(cbind(rowSums(Summer_Photic), rowSums(Fall_Photic),
                         rowSums(Summer_100m),rowSums(Fall_100m))>0)
names(comb)<-c("Summer Photic Zone", "Fall Photic Zone",
               "Summer 100 m","Fall 100 m")

venn(comb, borders=T, ggplot = T,col = "navyblue",zcolor = "bw",box = F)+
  geom_path(size=10)


ggsave("./Figures/Venn_Summer_Fall.pdf")
