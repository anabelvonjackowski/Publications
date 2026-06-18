############################################################################################
# MOLA RDA (Redundancy Analysis)
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Last Update: June 2024 from E. Fadeev et al. (2018)
############################################################################################

library("ggplot2"); packageVersion("ggplot2")   # For plotting the final multivariate ordination map
library("ggpmisc"); packageVersion("ggpmisc")   # For extended text/numerical plotting accents
library("phyloseq"); packageVersion("phyloseq") # For managing and subsetting target microbial data structures
library("vegan"); packageVersion("vegan")       # For ecological multivariate ordination transformations (e.g., rda)
library("cowplot"); packageVersion("cowplot")   # For streamlined theme layouts and canvas aggregations
library("dplyr"); packageVersion("dplyr")       # For structural matrix joining and pipes
library(DESeq2); packageVersion("DESeq2")       # For variance stabilizing transformations (vst) on count data
library(ggrepel)                                # For preventing overlapping label geometries

############################################################################################
## Preprocess bacterial OTU by prevalence of each taxa
############################################################################################
# Subsetting raw sequence pools via abundance constraints
pseq.abs = filter_taxa(all_data, function(x) sum(x > 3) > (0.03 * length(x)), TRUE)
pseq.abs.surface.subset <- subset_samples(pseq.abs, Depth_m == "5 m" & Filter.size == "FL")

# [Commented Out] Spatial and sample outlier masking for structural environmental cross-comparisons
# pseq.abs.surface.subset.subset <- subset_samples(pseq.abs, Filter.size == "0.2" & Sample_ID != "sa2178" & Sample_ID != "sa2174" &
#                                             Sample_ID != "sa2176" & Sample_ID != "sa2172" & 
#                                             Sample_ID != "sa2128" & Sample_ID != "sa2136" &
#                                             Sample_ID != "sa2110" & Sample_ID != "sa2108" & Sample_ID != "sa2104" & Sample_ID != "sa2158")
# pseq.abs.surface.subset@sam_data[,c(1:11, 88)]

# Clean target subset focusing strictly on 0.2 µm filter fractions at 5m surface depths
pseq.abs.surface.subset <- subset_samples(pseq.abs, 
                                          Depth_m == "5 m" & Filter.size == "0.2" &
                                            Sample_ID != "sa2176" & Sample_ID != "sa2172" & 
                                            Sample_ID != "sa2128" & Sample_ID != "sa2136" &
                                            Sample_ID != "sa2110")

# Determine global occurrence tracking (prevalence matrix) for each distinct taxon
prev0 = apply(X = otu_table(pseq.abs.surface.subset),
              MARGIN = ifelse(taxa_are_rows(pseq.abs.surface.subset), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

## Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq.abs.surface.subset)

## Execute prevalence filter, using `prune_taxa()` function
BAC_pruned <-  prune_taxa((prev0 > prevalenceThreshold), pseq.abs.surface.subset)

## Initiate Variance Stabilized Transformation (VST) framework via DESeq2 engine conversion
BAC_pruned.dds <- phyloseq_to_deseq2(BAC_pruned, ~1)
varianceStabilizingTransformation(BAC_pruned.dds, blind = TRUE, fitType = "parametric")
BAC_pruned.dds <- estimateSizeFactors(BAC_pruned.dds)

# Operational Error Bypass Layer: If log geometric mean fails due to zero-dominant matrices
dds <- BAC_pruned.dds[ rowSums(counts(BAC_pruned.dds)) > x, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
BAC_pruned.dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
varianceStabilizingTransformation(BAC_pruned.dds, blind = TRUE, fitType = "parametric")

# Extract final transformed variance matrix fields
BAC_pruned.dds <- estimateDispersions(BAC_pruned.dds)
otu.vst <- getVarianceStabilizedData(BAC_pruned.dds)

## QA Step: Confirm that the dimensions of the transformed matrix match the original pruned phyloseq object
dim(otu.vst)
dim(otu_table(BAC_pruned))

# Re-inject the variance stabilized data matrix back into a new companion phyloseq skeleton
BAC_pruned.vst<-BAC_pruned
otu_table(BAC_pruned.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

BAC_pruned.vst

## Assign workspace structures and filter unobserved operational taxonomic groups
BAC.no.na <- BAC_pruned.vst 

# BAC.no.na@sam_data$Depth <- BAC.no.na@sam_data$Depth
# %>% subset_samples(!is.na(Temperature_C))

# Remove taxa that display a zero sum tracking record across remaining sample objects
BAC.no.na <- prune_taxa(taxa_sums(BAC.no.na)>0,BAC.no.na)

############################################################################################
## Extract and Scale Environmental Parameters
############################################################################################
# [Commented Out] Selection matrix alternative 1: General Core physical metrics
# BAC.env <- data.frame(sample_data(BAC.no.na))[c("Temperature_C","Salinity_psu",
#                                                 "Chla_ugL","DOC_umolL")]

# [Commented Out] Selection matrix alternative 2: Free-living baseline amino acid arrays
# BAC.env <- data.frame(sample_data(BAC.no.na))[c("LAsp_nmolL","LThr_nmolL",
#                                                 "LArg_nmolL",
#                                                 "LHis_nmolL",
#                                                 "LLeu_nmolL",
#                                                 "DAsp_nmolL", "DAla_nmolL", "DLeu_nmolL")]

# Selection matrix alternative 3: Active environment parameters and pigment indicators
BAC.env <- data.frame(sample_data(BAC.no.na))[c("Temperature_C","Salinity_psu",
                                                "Chla_ugL",
                                                "DON_umolL",
                                                "DinoflagellatesA",
                                                "DAsp_nmolL", "DAla_nmolL")]

# Normalize environmental variable vectors to ensure uniform tracking scales across models
BAC.env <- as.data.frame(scale(BAC.env,center = FALSE, scale = TRUE))

## Extract transposed community data tracking tables for execution loops
BAC.otu <- t(otu_table(BAC.no.na))

# Compute Redundancy Analysis (RDA) formula model
BAC.rda.all <- rda(BAC.otu ~ ., data = BAC.env) 

############################################################################################
## Process Ordination Score Coordinates
############################################################################################
BAC.rda.scores <- vegan::scores(BAC.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC.rda.sites <- data.frame(BAC.rda.scores$sites)
BAC.rda.species <- data.frame(BAC.rda.scores$species*1.5)
BAC.rda.species.merged <- merge(BAC.rda.species,ampvis1, by = "row.names")

# Coordinate mapping keys for data frame table joining routines
BAC.rda.sites$Sample_ID <- as.character(rownames(BAC.rda.sites))
sample_data(BAC.no.na)$Sample_ID = as.character(rownames(sample_data(BAC.no.na)))
BAC.rda.sites <- BAC.rda.sites %>%
  left_join(sample_data(BAC.no.na))

# Format and isolate multi-variable directional vector arrows
BAC.rda.arrows<- BAC.rda.scores$biplot*5
colnames(BAC.rda.arrows)<-c("x","y")
BAC.rda.arrows <- as.data.frame(BAC.rda.arrows)

# Compute eigenvalues to map total variation percentages explained per ordination axis
BAC.rda.all$CCA$eig
BAC.rda.evals <- 100 * (BAC.rda.all$CCA$eig / sum(BAC.rda.all$CCA$eig))

# Set base graphics parameters
theme_set(theme_classic())
# theme_set(theme_bw())

# Configure categorical plotting vector levels
# BAC.rda.sites$Depth <- factor(BAC.rda.sites$Depth, levels=c("Surface","DCM","Below DCM", "100m"))
BAC.rda.sites$Month <- factor(BAC.rda.sites$Month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov"))

############################################################################################
## Generate and Export the Final ggplot2 Ordination Visual Layer
############################################################################################
ggplot() +
  # Map transformed tracking point geometries
  geom_point(data = BAC.rda.sites, aes(x = RDA1, y = RDA2, color = Depth_m), size = 4, shape = 16) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC.rda.evals[1], 2)),y = sprintf("RDA2 [%s%%]", round(BAC.rda.evals[2], 2))) +
  scale_color_manual(values = c("#F3DF6C","grey50"))+
  # scale_color_manual(values = c("Jan"="#3B9AB2", "Feb"="#0B775E","Mar"="#3db229","Apr"="#F3DF6C","May"="#DC863B","Jun"="#6d4401", "Jul"="#C7B19C",
  #                              "Aug"="#E6A0C4","Sep"="#F8AFA8","Nov"="grey50"))+
  # scale_x_reverse()+
  # scale_x_continuous(limits = c(-6,6))+
  # scale_y_continuous(trans = "reverse")+
  
  # Draw environmental variable directional load segments
  geom_segment(data=BAC.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  
  # Map environmental variable labels at a slight offset layout configuration multiplier
  geom_text(data=as.data.frame(BAC.rda.arrows*1.1), aes(x, y, label = rownames(BAC.rda.arrows)),color="black",alpha=1, fontface = "bold")+
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text= element_text(size = 12),
        plot.title = element_text(size=12, lineheight=.8, face="bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        plot.background = element_rect(color = "transparent"))

# Export high-resolution landscape publication graphic matching custom millimeter presets
