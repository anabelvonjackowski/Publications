##########################
# MOLA 16S Amplicon Sequencing Profile Processing
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Last Update: June 2024
##########################

library(phyloseq); packageVersion("phyloseq") # For managing and subsetting multi-table amplicon datasets
library(dplyr)                                 # For structural data pipelines and vector alignment
library("readxl")                              # For importing companion environmental sheet metadata

############################################################################################
###  PHYLOSEQ LOAD & PREPROCESSING  ###
############################################################################################
# Load Operational Taxonomic Unit (OTU) matrix
# OTU <- readRDS('./seqtab.rds')
# OTU <- read.csv("../MOLA_16S_files/OTU.csv",header = TRUE)

# Read raw OTU table and transpose to orient samples as rows and amplicon variants as columns
OTU <- t(read.csv("../MOLA_16S_files/OTU.csv", row.names = 1, header = TRUE))
# duplicated(OTU) | duplicated(OTU, fromLast = TRUE)

# write.table(t(as(otu_table(ps), "matrix")), "./OTU_table.txt", sep ="\t", quote = F)
# TAX <- readRDS('./taxa.rds')

# Load the corresponding taxonomic assignment matrix
TAX <- read.csv("../MOLA_16S_files/TAX.csv", row.names = 1, header = TRUE)
# write.table(as(tax_table(ps), "matrix"), "./TAX_table.txt", sep ="\t", quote = F)

# Import environmental metadata matrices, capturing rows 1-91 while stripping non-assigned entries
ENV <- subset(as.data.frame(read_excel('./mola_data.xlsx')[c(1:91),]),Sample_ID != "na")
rownames(ENV) <- ENV$Sample_ID
# ENV <- as.data.frame(read_excel('./molaenv.xlsx'))

# Build formal phyloseq-formatted OTU component tracking rows as taxa
otu <- otu_table(OTU, taxa_are_rows = TRUE)

# Quality Assurance Check: Ensure exact match symmetry between metadata sample keys and sequence tables
all.equal(sample_names(otu), ENV$Sample_ID)
# duplicated(OTU) | duplicated(OTU, fromLast = TRUE)

tax <- tax_table(TAX)
rownames(tax) <- rownames(otu)

# Construct comprehensive unified multi-table phyloseq container
all_data <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), 
                     sample_data(ENV),
                     tax_table(tax))

# Standardize classification schema header naming conventions
colnames(all_data@tax_table)<- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Remove specific aberrant or outlier samples
all_data <- subset_samples(all_data, Sample_ID != "S302FL"& Sample_ID != "S175FL")
# S302FL is strangely different from rest, , S175FL probably misplaced 3.0 µm filter????

# Apply taxonomic masking to omit eukaryotic contaminants (mitochondria and chloroplast sequences)
all_data1 <- subset_taxa(all_data, Family != "Mitochondria" & Order != "Chloroplast")
# View(all_data1@tax_table)
all_data <- all_data1

# Housekeeping: Purge unlinked workspace variables
rm(otu, tax)
# all_data_surface = subset_samples(all_data, Depth_m == "5 m")

############################################################################################
###  AMPVIS LOAD & SYNTAX INTEGRATION  ###
# Technical Documentation: https://kasperskytte.github.io/ampvis2/reference/index.html
############################################################################################
library(ampvis2)   # For stylized exploration of amplicon abundance profiles
library(lubridate) # For handling calendar string definitions (e.g., ymd)

# Low-abundance variance filter: Retain taxa presenting >3 counts in at least 3% of all samples
pseq.abs = filter_taxa(all_data, function(x) sum(x > 3) > (0.03 * length(x)), TRUE)

# Calculate internal sample proportions (relative abundance scaled to 100%)
pseq.rel = transform_sample_counts(pseq.abs, function(x) x / sum(x) * 100) 

# Coerce phyloseq data structures into structured flat dataframes for ampvis2 import compatibility
ampvis1 <- data.frame(OTU = rownames(phyloseq::otu_table(pseq.abs)@.Data),
                      phyloseq::otu_table(pseq.abs)@.Data,
                      phyloseq::tax_table(pseq.abs)@.Data,
                      check.names = FALSE)

# Extract and isolate metadata parameters
ampvis.env <- data.frame(phyloseq::sample_data(pseq.abs), check.names=F)
# ampvis.env <- cbind(Sample_ID = rownames(ampvis.env), ampvis.env)

# Bind structural inputs into dedicated native ampvis2 data objects
ampvis <- amp_load(ampvis1, ampvis.env)
ampvis$metadata$Year <- as.factor(ampvis$metadata$Year)
ampvis$metadata$Filter.size <- as.factor(ampvis$metadata$Filter.size)
ampvis$metadata$Month <- factor(ampvis$metadata$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

# Clear local context matrix memory footprint
rm(ampvis.env)

# Spatial and Seasonal Sample Subsetting:
ampvis_fL_subset <- amp_subset_samples(ampvis, Depth_m == "500 m" & Filter.size == "FL" & Season == "Spring")
# Spring FL 5m (n=3), Summer FL 5m (n=7), Fall FL 5m (n=3), Winter FL 5m (n=3)
# Spring FL 500m (n=3), Summer FL 5m (n=6), Fall FL 5m (n=4), Winter FL 5m (n=2)

library(readxl)
library(tidyverse)
library(lubridate)
library(anytime)
library(ampvis2)
library(phyloseq)

# Construct deep-sea environmental sample matrix targets
ampvis_fL <- amp_subset_samples(ampvis, Filter.size == "FL")
# ampvis_fL_surf <- amp_subset_samples(ampvis, Depth_m == "5 m")
ampvis_subset <- amp_subset_samples(ampvis, Depth_m == "500 m" & Filter.size == "FL")

############################################################################################
###  PLOTTING DATA LAYERS (HEATMAPS & TIMESERIES)
############################################################################################
# Generate Heatmap 1: Aggregated at Order classification level
amp_heatmap(ampvis_subset,
            tax_add = "Class",
            tax_aggregate = "Order",
            group_by = "JulianDate",
            # facet_by = "Year",
            tax_show = 30,
            plot_values = TRUE,
            showRemainingTaxa = TRUE,
            plot_colorscale = "log10")+
  theme( legend.position = "right")

# Generate Heatmap 2: Detailed perspective aggregated at Genus classification level
amp_heatmap(ampvis_subset,
            tax_add = "Class",
            tax_aggregate = "Genus",
            group_by = "JulianDate",
            # facet_by = "Year",
            tax_show = 100,
            plot_values = TRUE,
            showRemainingTaxa = TRUE,
            plot_colorscale = "log10")+
  theme( legend.position = "right")
# 
# Save Genus distribution blueprint profiles
ggsave("./Figures_16S/FL_Deep_Ampvis_Genus_100.png", height = 210*2, width = 350, dpi = 300, units = "mm") 

# Generate Abundance Timeseries Plots
amp_timeseries(ampvis_subset,
               time_variable = "JulianDate",
               # group_by = "Filter_size",
               tax_add = "Class",
               tax_aggregate = "Genus",
               # split = TRUE,
               scales = "free_y",
               tax_show = 10
)+
  # guides(color=guide_legend(title="Filter Size"))+
  # scale_color_manual(values = c("#046C9A", "#D69C4E"))+
  
  # Custom color assignment across dominant bacterial lines
  scale_color_manual(values = c("#455285",  "#AAB7E8","#1C5A62", "#3C6D56", "#687B3E", "#7A711F","#F8A17B", "#FDB7BC", 
                                "#E1AF00", "#DD8D29" ,"#E2D200","#F2AD00", "grey70"))+
  # scale_color_manual(values = c("#455285",  "#AAB7E8","#D3E0FA","#1C5A62", "#3C6D56", "#3C5600","#687B3E", "#7A711F","#9D892B", "#B79A5E","#D29343", "#F8A17B", "#FDB7BC", "#FACCFA","#F1CEA4", 
  # "#E1AF00", "#DD8D29" ,"#E2D200","#F2AD00", "grey70"))+
  geom_line(size=2)+
  geom_point(color = "black")+
  geom_vline(xintercept = 62)+  # Timeline marker 1
  geom_vline(xintercept = 112)+ # Timeline marker 2
  # scale_x_continuous(expand = c(0,0), limits = c(0,180))+
  # scale_x_date(expand = c(0,0),date_labels = "%Y",limits = as.Date(c('2019-01-01','2021-12-31')),
  # date_breaks = "1 year", date_minor_breaks = "2 months")+
  scale_y_continuous(expand = c(0,0), limits = c(0,60))+
  theme(axis.title=element_text(size=12),
        panel.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        panel.grid.major.x =  element_line(color = "grey50"),
        panel.grid.minor.x =  element_line(color = "white"))
# 
# Save community structural timeline charts
ggsave("./Figures_16S/FL_Surface_Timeseries_Genus.png", width = 10)

#####################################
# Barplot & Reshaping Matrices
#####################################
library(phyloseq)

# Melt ampvis structures into long formats optimized for downstream ggplot layouts
ampvis_long <- amp_export_long(
  ampvis_fL,
  metadata_vars = colnames(ampvis_fL$metadata),
  tax_levels = colnames(ampvis_fL$tax))

ampvis_long <- amp_export_long(ampvis_fL, metadata_vars = colnames(ampvis_fL$metadata),
                               tax_levels = c("OTU", "Class"))

pseq_subset <- subset_samples(pseq.abs, Filter.size == "FL") #& Year != "2020"
# pseq_subset <- subset_samples(pseq.abs, Filter.size == "3" & Depth_m == "5 m") #& Year != "2020"
# pseq_subset <- subset_samples(pseq.abs, Filter.size == "0.2" & Depth_m == "500 m")

# Agglomerate and condense tax profiles by Genus designations
genusabundance <- pseq.abs %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()  %>%                                          # Melt to long format
  filter(Abundance > 0.05) %>%                           # Filter out low abundance taxa
  arrange(Genus) 
genusabundance$Abundance <-genusabundance$Abundance*100

