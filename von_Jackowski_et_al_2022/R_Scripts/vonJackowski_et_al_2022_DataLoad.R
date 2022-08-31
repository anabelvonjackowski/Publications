
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## Data Load
############################################################################################

## This script: load and format amplicons (from dada2) and metadata

# Set working directory; save and load
# setwd("/")
# load("Pebcao.Rdata")

# Load packages
library(gtools)
library(phyloseq)
library(ampvis2)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(gtools)

# if you have problems:
library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
library(ampvis2); packageVersion("ampvis2") #devtools::install_github("MadsAlbertsen/ampvis2", force = TRUE)

############################################################################################
###  RAW DATA -- LOAD AND FORMATTING ###
############################################################################################

## ASV IMPORT and TAXONOMY ##

# Read ASVs
ASV <- read.table(
  "ASV.txt",
  h = T, 
  row.names=1,
  sep = "\t")

# Edit taxonomy in Excel: replace NA with last taxrank plus "unclassified" 
TAX <- read.table(
  "tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

############################################################################################

## METADATA ##

# Import metadata (season, nutrients etc)
ENV <- read.table("metadata.txt", h = T, sep = "\t", stringsAsFactors = F)
row.names(ENV) = ENV$Fram_ID
str(ENV)

# Generate factors for parameter of interest
ENV.WM <- ENV %>%
  na.omit(dplyr::select(Station,Season, Temperature_C,Salinity_PSU))
ENV$Depth <- factor(ENV$Depth, levels = c(
  "Surface", "DCM", "BDCM", "100m"))
ENV$Cruise <- factor(ENV$Cruise, levels = c(
  "PS114", "MSM77"))
ENV$Season <- factor(ENV$Season, levels = c(
  "Summer", "Fall"))

# rename BAC 
ASV <- ASV[,mixedsort(names(ASV))]
colnames(ASV) = ENV$Fram_ID


############################################################################################
###  PHYLOSEQ LOAD  ###
############################################################################################

asv = otu_table(ASV, taxa_are_rows = TRUE)
tax = tax_table(TAX)
rownames(tax) <- rownames(asv)

pseq <- phyloseq(otu_table(asv, taxa_are_rows = FALSE), 
                 sample_data(ENV), 
                 tax_table(tax))

# Fix tax-IDs
colnames(pseq@tax_table)<- c("Kingdom","Phylum","Class","Order","Family","Genus")

# remove temporary files
rm(asv, tax)

####################################

## Filter: only ASV with >3 counts in >3% of samples; transform to rel. ab
pseq.abs = filter_taxa(pseq, function(x) sum(x > 3) > (0.03 * length(x)), TRUE)
pseq.rel = transform_sample_counts(pseq.abs, function(x) x / sum(x) * 100) 

## Agglomerate taxranks
class.abs <- tax_glom(pseq.abs, taxrank = "Class")
class.rel <- tax_glom(pseq.rel, taxrank = "Class")

family.abs <- tax_glom(pseq.abs, taxrank = "Family")
family.rel <- tax_glom(pseq.rel, taxrank = "Family")

genus.abs <- tax_glom(pseq.abs, taxrank = "Genus")
genus.rel <- tax_glom(pseq.rel, taxrank = "Genus")

## Convert filtered abs/rel to dataframes
filt.rel = as(otu_table(pseq.rel), "matrix") 
filt.rel = as.data.frame(filt.rel)  
filt.abs = as(otu_table(pseq.abs), "matrix") 
filt.abs = as.data.frame(filt.abs) 
filt.tax = as(tax_table(pseq.rel), "matrix")  
filt.tax = as.data.frame(filt.tax)

# Merge and export
ASV.rel <- cbind(filt.rel, filt.tax)
ASV.abs <- cbind(filt.abs, filt.tax)
# write.table(ASV.rel, file="ASV.rel.txt", sep="\t", row.names = T)
# write.table(ASV.abs, file="ASV.abs.txt", sep="\t", row.names = T)

# Remova temporary datasets
rm(filt.rel, filt.abs, filt.tax)

############################################################################################
###  AMPVIS LOAD  ###
############################################################################################

# pseq.abs_subset <- subset_samples(pseq.abs, Season=="Summer")
# pseq.abs <-pseq.abs_subset

ampvis1 <- data.frame(OTU = rownames(phyloseq::otu_table(pseq.abs)@.Data),
   phyloseq::otu_table(pseq.abs)@.Data,
   phyloseq::tax_table(pseq.abs)@.Data,
   check.names = FALSE)

# Add "Species" column to BAC data (required >> copy "Genus" column)
ampvis1$Species <- ampvis1$Genus

# Extract metadata from phyloseq; format and combine
ampvis.env <- data.frame(phyloseq::sample_data(pseq.abs), check.names=F)
ampvis.env <- cbind(Sample = rownames(ampvis.env), ampvis.env) 
ampvis <- amp_load(ampvis1, ampvis.env)

                   
# remove temporary datasets
rm(ampvis.env)

#############################################################################################

# save.image("Pebcao.Rdata")

