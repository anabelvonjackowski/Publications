# © This script is modified from Eduard Fadeev @eddie1986
#####################################
#setup
#####################################
setwd("~/Documents/...")
rm(list = ls()) ##Clears entire workspace 
#Load data
load("./x.RData")
#####################################
#Load libraries
#####################################
library(phyloseq)
library(ggplot2)
library(DESeq2)
#####################################
#Plot total number of reads and OTUs per sample
#####################################
readsumsdf <- data.frame(nreads = sort(taxa_sums(all_data), TRUE), sorted = 1:ntaxa(all_data), 
                         type = "OTUs")
readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(all_data), 
                                                         TRUE), sorted = 1:nsamples(all_data), type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("./figures/readcounts_V4V5.pdf")
###########################
# OTU overview
###########################
OTU_over <- all_data@sam_data[,1:3]
OTU_over$Total.reads <- colSums(all_data@otu_table)       # number of Sequences per sample
OTU_over$Total.OTUs <- colSums(all_data@otu_table!=0)     # number of OTUs per sample
write.table(OTU_over,"./OTUs_counts_overview.tsv", sep=" ", row.names = F)
###########################
# OTU variation
###########################
hist(log10(apply(otu_table(all_data), 1, var)),
     xlab="log10(variance)", breaks=50,
     main="OTUs variance")
ggsave("./figures/OTUVarience_V4V5.pdf")
#####################################
#Rarefy the dataset by the smallest sample
#####################################
set.seed(12345)
all_data.rare <- rarefy_even_depth(all_data, sample.size = min(sample_sums(all_data)),
                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
all_data.rare.10K <- rarefy_even_depth(all_data, sample.size = 10000,
                                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#####################################
#Remove taxa not seen more than 3 times in at least 20% of the samples
#####################################

all_data.filt3 <- filter_taxa(all_data, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

#####################################
#Standardize abundances to the median sequencing depth
#####################################
total <- median(sample_sums(all_data))
standf <- function(x, t=total) round(t * (x / sum(x)))
all_data.median <- transform_sample_counts(all_data, standf)
################################################################################
# Variance stabilized Transformation
################################################################################

all_data.sub <- subset_samples(all_data.rare)

all.dds <- phyloseq_to_deseq2(all_data.sub, ~1)
all.dds <- estimateSizeFactors(all.dds)
all.dds <- estimateDispersions(all.dds, fitType = "local")
data.vst <- getVarianceStabilizedData(all.dds)

#make sure that the dimentions aof the OTU table and the DEseq object are matching
dim(data.vst)
dim(all_data.sub@otu_table)

all_data.vst<-all_data.sub
otu_table(all_data.vst)<- otu_table(data.vst, taxa_are_rows = TRUE)

################################################################################
# Regularized Log transformation (rLog) Transformation
################################################################################
#rlog transformation
all.dds <- phyloseq_to_deseq2(all_data.sub, ~1)
data.rld <- rlog(all.dds,  blind=TRUE, fitType = "local")
rlogMat <- assay(data.rld)
rownames(rlogMat)<-rownames(data.rld)
colnames(rlogMat)<-colnames(data.rld)

all_data.rld<-all_data.sub
otu_table(all_data.rld)<- otu_table(rlogMat, taxa_are_rows = TRUE)

################################################################################
# Centered Log-Ratio (CLR) Transformation
################################################################################
gm_mean = function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
clr = function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

all_data.clr <- transform_sample_counts(all_data.sub, fun = clr)
#save 
save.image(file="../Preprossessing_data.RData")
