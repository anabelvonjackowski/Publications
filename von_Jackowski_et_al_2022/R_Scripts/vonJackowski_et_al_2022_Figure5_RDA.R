
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## #Forward selection of explanatory variables
############################################################################################
# based on the script by Eddie Fadeev @edfadeev

library("ggplot2"); packageVersion("ggplot2")
library("ggpmisc"); packageVersion("ggpmisc")
library("phyloseq"); packageVersion("phyloseq")
library("vegan"); packageVersion("vegan")
library("cowplot"); packageVersion("cowplot")
library("dplyr"); packageVersion("dplyr")
library(DESeq2)
library(ggrepel)

############################################################################################
## Preprocess bacterial OTU by prevalence of each taxa
############################################################################################

prev0 = apply(X = otu_table(pseq.abs),
              MARGIN = ifelse(taxa_are_rows(pseq.abs), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

## Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq.abs)

## Execute prevalence filter, using `prune_taxa()` function
BAC_pruned <-  prune_taxa((prev0 > prevalenceThreshold), pseq.abs)

## Variance stabilized Transformation after prevalance
BAC_pruned.dds <- phyloseq_to_deseq2(BAC_pruned, ~1)
varianceStabilizingTransformation(BAC_pruned.dds, blind = TRUE, fitType = "parametric")
BAC_pruned.dds <- estimateSizeFactors(BAC_pruned.dds)
BAC_pruned.dds <- estimateDispersions(BAC_pruned.dds)
otu.vst <- getVarianceStabilizedData(BAC_pruned.dds)

## make sure that the dimentions aof the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(BAC_pruned))

BAC_pruned.vst<-BAC_pruned
otu_table(BAC_pruned.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

BAC_pruned.vst

## subset by fraction and remove NAs in metadata
BAC.no.na <- BAC_pruned.vst 

# BAC.no.na@sam_data$Depth <- BAC.no.na@sam_data$Depth
# %>% subset_samples(!is.na(Temperature_C))

#remove unobserved OTU
BAC.no.na <- prune_taxa(taxa_sums(BAC.no.na)>0,BAC.no.na)

############################################################################################
## Extract the environmental parameters
############################################################################################

BAC.env <- data.frame(sample_data(BAC.no.na))[c("Galactose","Fucose","Rhamnose","Tyrosine","Glycine",
                                                      "Glucose",	"Aspartic.Acid","Arabinose",	"Alanine","Galactosamine","Galacturonic.Acid","Mannose_Xylose",	"Glucuronic.Acid" ,"Glucosamine","Threonine", "Arginine",
                                                      "Glutamic.Acid","GABA","Isoleucine",	"Phenylalanine",	"Valine",		"Serine",	"Leucine")]
BAC.env <- as.data.frame(scale(BAC.env,center = FALSE, scale = TRUE))

## Extract OTU tables from Phyloseq object
BAC.otu <- t(otu_table(BAC.no.na))

BAC.rda.all <- rda(BAC.otu ~ ., data = BAC.env) # model including all variables 

############################################################################################
## Generate an RDA plot 
############################################################################################
BAC.rda.scores <- vegan::scores(BAC.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC.rda.sites <- data.frame(BAC.rda.scores$sites)
BAC.rda.species <- data.frame(BAC.rda.scores$species*1.5)
BAC.rda.species.merged <- merge(BAC.rda.species,ampvis1, by = "row.names")

BAC.rda.sites$Sample_ID <- as.character(rownames(BAC.rda.sites))
sample_data(BAC.no.na)$Sample_ID <- as.character(rownames(sample_data(BAC.no.na)))
BAC.rda.sites <- BAC.rda.sites %>%
  left_join(sample_data(BAC.no.na))

#Draw biplots
BAC.rda.arrows<- BAC.rda.scores$biplot*5
colnames(BAC.rda.arrows)<-c("x","y")
BAC.rda.arrows <- as.data.frame(BAC.rda.arrows)
BAC.rda.evals <- 100 * (BAC.rda.all$CCA$eig / sum(BAC.rda.all$CCA$eig))

#set plots theme
theme_set(theme_classic())
# theme_set(theme_bw())

BAC.rda.sites$Depth <- factor(BAC.rda.sites$Depth, levels=c("Surface","DCM","Below DCM", "100m"))
BAC.rda.sites$Station <- factor(BAC.rda.sites$Station, levels=c("N5","N4","N3", "S3","HG4","HG1/2","SV4","SV2"))

ggplot() +
  geom_point(data = BAC.rda.sites, aes(x = RDA1, y = RDA2, color = Depth), size = 4, shape = 16) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC.rda.evals[1], 2)),y = sprintf("RDA2 [%s%%]", round(BAC.rda.evals[2], 2))) +
  scale_color_manual(values = c("#b34f96","#71ae5d","#6a65bd","#bc4949"))+
  # scale_x_reverse()+
  scale_x_continuous(limits = c(-6,6))+
  scale_y_continuous(trans = "reverse")+
  geom_segment(data=BAC.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC.rda.arrows*1.1), aes(x, y, label = rownames(BAC.rda.arrows)),color="black",alpha=1, fontface = "bold")+
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text= element_text(size = 12),
        plot.title = element_text(size=12, lineheight=.8, face="bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1), 
        plot.background = element_rect(color = "transparent"))

ggsave("./Figures/Betadiversity_RDA_ASV.pdf", width = 9, height = 7)
