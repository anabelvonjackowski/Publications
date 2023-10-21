# © This script is modified by Anabel von Jackowski @avonjackowski
# based on https://github.com/edfadeev/Bact-comm-PS85
#####################################
#load libraries
#####################################
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("phyloseq"); packageVersion("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("gage"); packageVersion("gage")
library("dplyr"); packageVersion("dplyr")

#####################################
# Set Basic Parameters
#####################################
#set plots theme
theme_set(theme_classic())

#load dataset
BAC <- all_data
# BAC <- all_data_MUC

# prune dataset
BAC_prune <- prune_taxa(taxa_sums(BAC)>0,BAC)
# BAC_subset <- subset_samples(BAC_prune, Depth_s == "0-10m"  | Depth_s == "50 m")
BAC_subset <- subset_samples(BAC_prune, Depth_s == "50 m" | Depth_s == "Intermediate")
# BAC_subset <- subset_samples(BAC_prune, Depth_s == "Intermediate" | Depth_s == "Deep")
BAC_subset@sam_data$Depth_s
# BAC_subset <- subset_samples(BAC_prune, Layer != "Deep Layer (14-16cm)")
# BAC_subset <- subset_samples(BAC_prune, Layer == "Deep Layer (14-16cm)")
# BAC_subset@sam_data$Layer_g
# BAC_subset@sam_data$Location_NMDS


#define alpha
alpha_val <- 0.05
q_val <- 0.05

#####################################
##run DEseq2
####################################
BAC.ddsMat <- phyloseq_to_deseq2(BAC_subset, ~Depth_s)
# BAC.ddsMat <- phyloseq_to_deseq2(BAC_subset, ~Location_NMDS)
varianceStabilizingTransformation(BAC.ddsMat, blind = TRUE, fitType = "parametric")
BAC.ddsMat <- estimateSizeFactors(BAC.ddsMat)
BAC.ddsMat <- estimateDispersions(BAC.ddsMat)

BAC.DEseq <- DESeq(BAC.ddsMat, fitType="parametric")
BAC.DEseq.res <- results(BAC.DEseq)

#####################################
#create taxonomy db
#####################################
# Identify enriched taxa
#create taxonomy db
# tax_subset <- subset_samples(all_data, Depth_s == "0-10m"  | Depth_s == "50 m")
tax_subset <- subset_samples(all_data,   Depth_s == "Intermediate"  | Depth_s == "50 m")
# tax_subset <- subset_samples(all_data,   Depth_s == "Intermediate" | Depth_s == "Deep")
tax_subset@sam_data$Depth_s
# tax_subset <- subset_samples(all_data, Layer != "Deep Layer (14-16cm)")
# tax_subset <- subset_samples(all_data, Layer == "Deep Layer (14-16cm)")
# tax_subset@sam_data$Layer

BAC_tax <- as.data.frame(tax_table(tax_subset))
BAC_tax$OTU <- rownames(tax_table(tax_subset))
BAC_tax <- BAC_tax %>% mutate_if(is.factor, as.character)

# Extract the desired taxonomic level
# Genera <- unique(BAC_tax$Class)
# colnames(BAC_tax)[3] <- "id"

# Family <- unique(BAC_tax$Family)
# colnames(BAC_tax)[5] <- "id"

Genera <- unique(BAC_tax$Genus)
colnames(BAC_tax)[6] <- "id"

#generate list of OTU for each Taxa
OTU.gs <- list()

for (s in 1:length(Genera)){
  n <- Genera[s]
  Order_OTU <- subset(BAC_tax, id == n)
  OTU.gs[[n]] <- Order_OTU$OTU
  
}

#####################################
#Identify enriched families / genera
#####################################
deseq2.fc <- BAC.DEseq.res$log2FoldChange
names(deseq2.fc) <- rownames(BAC.DEseq.res)
exp.fc=deseq2.fc

fc.Order.p <- gage(exp.fc, gsets = OTU.gs, same.dir=TRUE, ref = NULL, samp = NULL)

#plot
enrch <- rbind(fc.Order.p$greater,fc.Order.p$less)
enrch <-  data.frame(enrch[enrch[,"q.val"]<0.05 &
                                   !is.na(enrch[,"q.val"]),])
enrch$id <- rownames(enrch)
# enrch <- unique(merge(enrch,BAC_tax[,c("id","Phylum")]))
enrch <- unique(merge(enrch,BAC_tax[,c("id","Class")]))
enrich <- enrch
#####################################
# Correct taxa (optional)
#####################################
# 0 to 50  Class #####
enrich$id <- gsub("BD2-11 terrestrial group", "Other",enrich$id)
enrich$id <- gsub("Clostridia", "Other",enrich$id)
enrich$id <- gsub("Kiritimatiellae", "Other",enrich$id)
enrich$id <- gsub("PAUC34f_unclassified", "Other",enrich$id)
enrich$id <- gsub("Pla3 lineage", "Other",enrich$id)
enrich$id <- gsub("Subgroup 21", "Other",enrich$id)
enrich$id <- gsub("Subgroup 6", "Other",enrich$id)
enrich$id <- gsub("Subgroup 5", "Other",enrich$id)
enrich$id <- gsub("TK10", "Other",enrich$id)
enrich$id <- gsub("TK17", "Other",enrich$id)
enrich$id <- gsub("Woesearchaeia", "Other",enrich$id)
enrich$id <- gsub("Actinobacteria", "Other",enrich$id)
enrich$id <- gsub("Margulisbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("OM190", "Other",enrich$id)
enrich$id <- gsub("P9X2b3D02", "Other",enrich$id)
enrich$id <- gsub("Acidobacteriia", "Other",enrich$id)
enrich$id <- gsub("Dehalococcoidia", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedentia", "Other",enrich$id)
enrich$id <- gsub("JG30-KF-CM66", "Other",enrich$id)
enrich$id <- gsub("Nitrospira", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedentia", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedentia", "Other",enrich$id)
enrich <- subset(enrich, id != "Other")

# 0-50 genus 
enrich$Class <- gsub("Subgroup 6", "Acidobacteria",enrich$Class)
enrich$id <- gsub("Actinomarinales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Alphaproteobacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Anaerolineaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Endoecteinascidia", "Other",enrich$id)
enrich$id <- gsub("Coxiella", "Other",enrich$id)
enrich$id <- gsub("EPR3968-O8a-Bc78_unclassified", "Other",enrich$id)
enrich$id <- gsub("Flavobacteriaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("FS140-16B-02 marine group", "Other",enrich$id)
enrich$id <- gsub("Ga0077536_unclassified", "Other",enrich$id)
enrich$id <- gsub("Gimesiaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedensaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Margulisbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("MB11C04 marine group", "Other",enrich$id)
enrich$id <- gsub("NB1-j_unclassified", "Other",enrich$id)
enrich$id <- gsub("Oligoflexaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("OM27 clade", "Other",enrich$id)
enrich$id <- gsub("P3OB-42_unclassified", "Other",enrich$id)
enrich$id <- gsub("P9X2b3D02_unclassified", "Other",enrich$id)
enrich$id <- gsub("PAUC26f", "Other",enrich$id)
enrich$id <- gsub("PAUC34f_unclassified", "Other",enrich$id)
enrich$id <- gsub("PB19_unclassified", "Other",enrich$id)
enrich$id <- gsub("Rubinisphaeraceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Sneathiellaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 5_unclassified", "Other",enrich$id)
enrich$id <- gsub("UBA10353 marine group_unclassified", "Other",enrich$id)
enrich$id <- gsub("WCHB1-41_unclassified", "Other",enrich$id)
enrich$id <- gsub("AEGEAN-169 marine group_unclassified", "Other",enrich$id)
enrich$id <- gsub("Woeseia", "Other",enrich$id)
enrich <- subset(enrich, id != "Other")

# 200 to 400  #####
# Class #####
enrich$id <- gsub("Pla3 lineage", "Other",enrich$id)
enrich$id <- gsub("TK17", "Other",enrich$id)

# Genus
enrich$Class <- gsub("Subgroup 6", "Acidobacteria",enrich$Class)
enrich$id <- gsub("Mitochondria_unclassified", "Other",enrich$id)
enrich$id <- gsub("OM43 clade", "Other",enrich$id)
enrich$id <- gsub("Planctomycetales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Planktomarina", "Other",enrich$id)
enrich$id <- gsub("SAR116 clade_unclassified", "Other",enrich$id)
enrich <- subset(enrich, id != "Other")

#400 to 600 ##### 
enrich$id <- gsub("Margulisbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 21", "Other",enrich$id)
enrich <- subset(enrich, id != "Other")

enrich$id <- gsub("Margulisbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Oligoflexaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("SAR11 clade_unclassified", "Other",enrich$id)
enrich <- subset(enrich, id != "Other")

# 
# Sediment Surf-Sub #####
enrich$id <- gsub("B2M28_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Kuenenbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Peregrinibacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Uhrbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Woykebacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Corynebacterium 1", "Other",enrich$id)
enrich$id <- gsub("Gracilibacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedensaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Leptospirillum", "Other",enrich$id)
enrich$id <- gsub("Rhodopirellula", "Other",enrich$id)
enrich$id <- gsub("Rubritalea", "Other",enrich$id)
enrich$id <- gsub("Sva0996 marine group", "Other",enrich$id)
enrich$id <- gsub("Thiotrichaceae_unclassified", "Other",enrich$id)
# enrich$id <- gsub("Rhodobacteraceae_unclassified", "Rarinimicrobia (SAR406 clade)_unclassified",enrich$id)
enrich <- subset(enrich, id != "Other")
# 
# Sediment Sub-Deep #########
# enrich$Class <- gsub("Subgroup 17", "Acidobacteria",enrich$Class)
enrich$id <- gsub("Subgroup 10", "Other",enrich$id)
enrich$id <- gsub("Cyclobacteriaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 6_unclassified", "Other",enrich$id)
enrich$id <- gsub("Dadabacteriales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Unknown Family_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Yanofskybacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Flavobacteriaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 9_unclassified", "Other",enrich$id)
enrich$id <- gsub("PAUC43f marine benthic group_unclassified", "Other",enrich$id)
enrich$id <- gsub("OM190_unclassified", "Other",enrich$id)
enrich$id <- gsub("Omnitrophicaeota_unclassified", "Other",enrich$id)
enrich$id <- gsub("AqS1", "Other",enrich$id)
enrich$id <- gsub("Arenicellaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Rhodobacteraceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Nitrosomonas", "Other",enrich$id)
enrich$id <- gsub("Anaerolineaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 2_unclassified", "Other",enrich$id)
enrich$id <- gsub("67-14_unclassified", "Other",enrich$id)
enrich$id <- gsub("Arcobacter", "Other",enrich$id)
enrich$id <- gsub("BRC1_unclassified", "Other",enrich$id)
enrich$id <- gsub("BD7-8_unclassified", "Other",enrich$id)
enrich$id <- gsub("Blastopirellula", "Other",enrich$id)
enrich$id <- gsub("Candidatus Berkiella", "Other",enrich$id)
enrich$id <- gsub("Candidatus Buchananbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Kerfeldbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Blastopirellula", "Other",enrich$id)
enrich$id <- gsub("Blfdi19_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Kuenenbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Peregrinibacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Woesebacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Candidatus Woykebacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Cellvibrionales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Colwelliaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("EC3_unclassified", "Other",enrich$id)
enrich$id <- gsub("Ekhidna", "Other",enrich$id)
enrich$id <- gsub("Halioglobus", "Other",enrich$id)
enrich$id <- gsub("BD1-7 clade", "Other",enrich$id)
enrich$id <- gsub("B2M28_unclassified", "Other",enrich$id)
enrich$id <- gsub("AT-s2-59_unclassified", "Other",enrich$id)
enrich$id <- gsub("AKAU3564 sediment group_unclassified", "Other",enrich$id)
enrich$id <- gsub("EV818SWSAP88_unclassified", "Other",enrich$id)
enrich$id <- gsub("Fibrobacteraceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Hydrogenedensaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Hyphomicrobiaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("IMCC26256_unclassified", "Other",enrich$id)
enrich$id <- gsub("JGI 0000069-P22_unclassified", "Other",enrich$id)
enrich$id <- gsub("Latescibacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Lentimicrobiaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Lentisphaera", "Other",enrich$id)
enrich$id <- gsub("Leptospirillum", "Other",enrich$id)
enrich$id <- gsub("Margulisbacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Marinoscillum", "Other",enrich$id)
enrich$id <- gsub("MB-A2-108_unclassified", "Other",enrich$id)
enrich$id <- gsub("Methylophagaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("MVP-88_unclassified", "Other",enrich$id)
enrich$id <- gsub("Mycobacterium", "Other",enrich$id)
enrich$id <- gsub("Nitrosococcaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Nitrospina", "Other",enrich$id)
enrich$id <- gsub("Nitrospinaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("NS9 marine group_unclassified", "Other",enrich$id)
enrich$id <- gsub("Oleiphilus", "Other",enrich$id)
enrich$id <- gsub("OM182 clade_unclassified", "Other",enrich$id)
enrich$id <- gsub("OM60(NOR5) clade", "Other",enrich$id)
enrich$id <- gsub("Paenibacillus", "Other",enrich$id)
enrich$id <- gsub("Parcubacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Poribacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Proteobacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Pseudoalteromonadaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Pseudomonas", "Other",enrich$id)
enrich$id <- gsub("Psychromonas", "Other",enrich$id)
enrich$id <- gsub("Rhizobiales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Rhodopirellula", "Other",enrich$id)
enrich$id <- gsub("Rickettsiaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Rickettsiales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Saccharimonadales_unclassified", "Other",enrich$id)
enrich$id <- gsub("Sandaracinaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Saprospiraceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("SG8-4_unclassified", "Other",enrich$id)
enrich$id <- gsub("SM1A02", "Other",enrich$id)
enrich$id <- gsub("SM23-31", "Other",enrich$id)
enrich$id <- gsub("Solirubrobacter", "Other",enrich$id)
enrich$id <- gsub("Spongiibacteraceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("SS1-B-04-55_unclassified", "Other",enrich$id)
enrich$id <- gsub("Subgroup 12_unclassified", "Other",enrich$id)
enrich$id <- gsub("Thiotrichaceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("UASB-TL25_unclassified", "Other",enrich$id)
enrich$id <- gsub("UBA10353 marine group_unclassified", "Other",enrich$id)
enrich$id <- gsub("wb1-A12", "Other",enrich$id)
enrich$id <- gsub("WCHB1-41_unclassified", "Other",enrich$id)
enrich$id <- gsub("Zixibacteria_unclassified", "Other",enrich$id)
enrich$id <- gsub("Otherceae_unclassified", "Other",enrich$id)
enrich$id <- gsub("Schekmanbacteria_unclassified", "Sarinimicrobia (SAR406 clade)_unclassified",enrich$id)
enrich <- subset(enrich, id != "OM60(NOR5) clade")
enrich <- subset(enrich, id != "Other")
#####################################
# Plot
#####################################
# enrch.plot <- 
  ggplot(data=enrich, aes(y=stat.mean , x=id, label = set.size))+ 
  scale_x_discrete(position="top", limits = rev, " ")+
  # ggtitle("Seamount versus Reference/Vent")+
  ylab("Fold Change (log2)")+ 
  # geom_point(size = 5, aes(colour = Phylum))+
  geom_point(size = 5, aes(colour = Class))+
  geom_text(size = 4, aes(y=stat.mean , x=id), nudge_y= 2.5, nudge_x= 0)+
  ylim(15,-15)+
  # ylim(-15,15)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = my_color_class) +
  theme(axis.title.x = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(angle = 90,size = 12,  hjust=0.95,vjust=0.5, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, vjust = 0.3,hjust = 0, color = "black"),
      # axis.text.x = element_text(angle = 90,size = 12, hjust=0,vjust=0.5), #if switch "x"
      strip.text = element_text(size = 12, angle = 90),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      # legend.position= "right",
      legend.position= "none",
      panel.background = element_rect(fill = 'white', colour = 'white'),
      panel.border = element_rect(color = "black", fill = NA, size = 1), 
      strip.background = element_rect(color = "black", size = 1),
      strip.text.y = element_text(size = 12))+
  # facet_wrap(~Depth)+
  coord_flip()

#####################################
# Save
#####################################
write.csv(enrich, "./Figure-7-enriched_genus_intermediate-deep.csv")
ggsave("./figures/Figure-7-enriched_genus_intermediate-deep.pdf", dpi = 300, device= "pdf",
       # height = 146,
       # height = 80,
       height = 63,
       width = 150,
       # width = 250,
       units = "mm")
