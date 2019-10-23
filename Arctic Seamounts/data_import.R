# © This script is modified from Eduard Fadeev @eddie1986
#####################################
#setup
#####################################
# setwd("~/Documents/...") #set working directory
rm(list = ls()) ##Clears entire workspace
sessionInfo()
#####################################
#Load libraries
#####################################
library(phyloseq)
source("../scripts/ReadAmplicon.R")
#####################################
#Parse dataset for Phyloseq
#####################################
# This is only for Bacteria

#!BE CAREFUL WITH YOUR INPUT!!!!!!
# X_B ... domain = "Bacteria"
# X_A ... domain = "Archaea"

#WATER v4v5
X_B <- ReadAmplicon(otu = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD_v4v5/OTU_contingency_table.csv", 
                    tax = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD_v4v5/amplicons_seeds_taxonomy.txt", 
                    silva = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/scripts/silva132_tax_ssu_curated.tsv", 
                    domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = T)
X_A <- ReadAmplicon(otu = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD_v4v5/OTU_contingency_table.csv", 
                    tax = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD_v4v5/amplicons_seeds_taxonomy.txt", 
                    silva = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/scripts/silva132_tax_ssu_curated.tsv", 
                    domain = "Archaea", silva.sep = "\t", singletons = F, unclassified = T, write.files = T)
#v3v4
X_B <- ReadAmplicon(otu = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD-sed_v3v4/OTU_contingency_table.csv", 
                    tax = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/CTD-sed_v3v4/amplicons_seeds_taxonomy.txt", 
                    silva = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/scripts/silva132_tax_ssu_curated.tsv", 
                    domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = T)
#v3v5 ARCHAEA
X_A <- ReadAmplicon(otu = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/Sed_v3v5/OTU_contingency_table.csv", 
                    tax = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/16S_pipeline/Silva132/Sed_v3v5/amplicons_seeds_taxonomy.txt", 
                    silva = "~/Documents/Lebenslauf/5_MarMic M.Sc/Master_Thesis/16S/scripts/silva132_tax_ssu_curated.tsv", 
                    domain = "Archaea", silva.sep = "\t", singletons = F, unclassified = T, write.files = T)

#Bacterial OTUs
OTU <- data.frame(X_B$OTU)
OTU <- data.frame(X_A$OTU)
# Do if you are analyzing the full microbial community - Bacteria & Archaea
# OTU_B <- data.frame(X_B$OTU)
# OTU_A <- data.frame(X_A$OTU)
# OTU <- rbind(OTU_A, OTU_B)

#Bacterial Taxonomy
TAX <- X_B$TAX
TAX <- X_A$TAX
## Do if you are analyzing the full microbial community - Bacteria & Archaea
# TAX_B <- X_B$TAX
# TAX_A<- X_A$TAX
# TAX <- rbind(TAX_A, TAX_B)

# Contextual data
#V3v4
ENV <- read.csv("./Metadata_V3V4.csv", sep = "," , h = T, row.names = 1, fill = T)
#V4V5
ENV <- read.csv("./Watercolumn_Metadata_V4V5.csv", sep = "," , h = T, row.names = 1, fill = T)
#V3V5
ENV <- read.csv("./Sediment_Metadata_V3V5.csv", sep = "," , h = T, row.names = 1, fill = T)

#####################################
# Check order of samples
#####################################
all.equal(rownames(OTU), rownames(TAX))

#####################################
#export the different tables for further use - have everything in the tables
#####################################
write.table(OTU, "OTU_table.txt", sep ="\t", quote = F)
write.csv(OTU, "OTU_table.csv")
write.table(TAX, "TAX_table.txt", sep ="\t", quote = F)
write.csv(TAX, "TAX_table.csv")
write.table(ENV, "ENV_table.txt", sep ="\t", quote = F)
# otu_tax <- merge(OTU, TAX)

#####################################
#creating Phyloseq dataset
#####################################
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV) #, c(as.factor(ENV$Depth),as.factor(ENV$Depth_m), as.factor(ENV$Depth_s)))
all.data <- phyloseq(OTU, TAX, meta)

#####################################
#Remove unobserved OTUs
#####################################
#Are there any OTUs included in this dataset that have no counted reads at all prior to preprocessing?
# any(taxa_sums(all.data) == 0)
# #How many?
# sum(taxa_sums(all.data) == 0)
# #remove unobserved OTUs
# all.data <- prune_taxa(taxa_sums(all.data)>0, all.data)
# #Verification of the removal (needs to be FALSE)
# any(taxa_sums(all.data) == 0)
# #####################################
#Fix taxonomy 
#####################################
colnames(all.data@tax_table)<- c("Domain","Phylum","Class","Order","Family","Genus")
all.data@tax_table[is.na(all.data@tax_table)]<- "Unclassified"

# all.data <- subset_samples(all.data, Device !="CTD")
# all.data2 <- subset_samples(all.data, Device =="MUC")

#Remove OTUs
#v3v4
badOTU = subset_samples(all.data, Station == 'negContr')
badOTUs = c("otu791", "otu747","otu2271","otu2934","otu2171","otu620","otu3181","otu4602","otu2445","otu10548","otu11761",
            "otu1369","otu1450","otu77","otu101","otu206","otu251","otu501","otu807","otu2295","otu5097","otu71628")
#v3v5
badOTU = subset_samples(all.data, Station == 'neg')
badOTUs = c("otu1", "otu2", "otu5", "otu8", "otu11", "otu22", "otu28","otu379")

#RUN to exclude negative control
allTaxa = taxa_names(all.data)
allTaxa <- allTaxa[!(allTaxa %in% badOTUs)]
ex1 = prune_taxa(allTaxa, all.data)

#v3v4
all_data <- subset_samples(ex1, Device =="MUC") 
#v3v5
all_data = subset_samples(ex1, Station != 'neg')
subset_taxa(all_data, is.na(Class))


# 
# #Subset by devide
# all_data_CTD <- subset_samples(goodOTUs, Device=="CTD")
# CTD.otu <- all_data_CTD@otu_table
# plot_bar(CTD.otu)
# write.csv(all_data_CTD@otu_table, "OTU_table_CTD.csv")
# write.csv(all_data_CTD@tax_table, "TAX_table_CTD.csv")

write.csv(all_data@otu_table, "OTU_table_MUC_v3v5.csv")
write.csv(all_data_MUC@tax_table, "TAX_table_MUC.csv")

#save image
save.image(file="./phyloseq_v3v4_raw.RData")