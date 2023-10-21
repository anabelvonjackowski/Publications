# © This script is written by Anabel von Jackowski @avonjackowski
#####################################
#Load
#####################################
library(phyloseq)
library(ampvis2)
library(dplyr)
library(DESeq2)
load("./data_import_v3v4.RData")
load("./data_import_v4v5.RData")
load("./data_import_v3v5.RData")

############################################################################################
###  AMPVIS LOAD  ###
############################################################################################
pseq.abs = all_data
# pseq.abs = filter_taxa(all_data, function(x) sum(x > 3) > (0.03 * length(x)), TRUE)
# pseq.rel = transform_sample_counts(pseq.abs, function(x) x / sum(x) * 100) 

ampvis1 <- data.frame(OTU = rownames(phyloseq::otu_table(pseq.abs)@.Data),
                      phyloseq::otu_table(pseq.abs)@.Data,
                      phyloseq::tax_table(pseq.abs)@.Data,
                      check.names = FALSE)

# Extract metadata from phyloseq; format and combine
ampvis.env <- data.frame(phyloseq::sample_data(pseq.abs), check.names=F)
ampvis.env <- cbind(Sample = rownames(ampvis.env), ampvis.env) 
ampvis <- amp_load(ampvis1, ampvis.env)


# remove temporary datasets
rm(ampvis.env)

############################################################################################
###  Plot Top Taxa using AMPVIS Package
############################################################################################
#v4v5 Water Column
# ampvis.subset <- amp_subset_samples(ampvis, Sample != "X27" & Sample != "X26" & Sample != "X25")
# ampvis.subset$metadata$Depth <- factor(ampvis.subset$metadata$Depth, levels = c("Surface", "200 m","400 m","600 m", "Bottom Water"))
ampvis$metadata$Depth_s <- factor(ampvis$metadata$Depth_s, levels = c("0-10m","50 m","Intermediate","Deep"))
# ampvis.subset$metadata$Location <- factor(ampvis.subset$metadata$Location, levels = c("Reference Site","Northern Seamount", "Central Seamount","Karasik Seamount"))

#Sediment
# ampvis.subset$metadata$Location <- factor(ampvis.subset$metadata$Location, levels = c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle","Karasik Seamount", "Southern Slope","Reference Site"))
# ampvis$metadata$Layer_g <- factor(ampvis$metadata$Layer_g, levels = c("SurfaceSubsurface","Deep"))
ampvis$metadata$Layer <- factor(ampvis$metadata$Layer, levels = c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))

amp_heatmap(data = ampvis,
            # tax_aggregate = "Class", showRemainingTaxa = TRUE,
            tax_aggregate = "Genus", showRemainingTaxa = TRUE,
            # tax_add = "Class",
            group_by = "Layer",
            # facet_by = "Layer_g",
            # group_by = "Location_NMDS",
            # group_by = "Depth_s",
            tax_show = 5,
            min_abundance = 1.0)+
  scale_fill_gradientn(colours = c("white", "#4363d8"))+
  theme(legend.position = "right")
ggsave("./figures/Figure-S7-enriched_class_abundance_v3v5.pdf", dpi = 300, device= "pdf", height = 70, width = 150, units = "mm")
# ggsave("./figures/Figure-S7-enriched_genus_abundance_v3v4.pdf", dpi = 300, device= "pdf", height = 150, width = 235, units = "mm")

## Order Genera by Class
# View(ampvis.subset[["tax"]])
amp_heatmap(data = ampvis.subset,
            # tax_add = "Phylum",
            tax_aggregate = "Genus", showRemainingTaxa = TRUE,
            tax_show = 5,normalise = T,plot_values = T,
            order_y_by = c("Woesearchaeia_unclassified",
                           "Thaumarchaeota_unclassified","Nitrosopumilaceae_unclassified","Candidatus Nitrosopumilus","Marine Benthic Group A_unclassified"),
            min_abundance = 0.1,plot_colorscale = "log10",plot_values_size = 4,
            # group_by = "Location",facet_by = "Depth",
            group_by = "Location",facet_by = "Layer",
            color_vector = c("grey","whitesmoke","white"))+
  scale_y_discrete(position = "right")+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        # panel.grid.major = element_line(color = "black"), # get rid of major grid
        # panel.grid.minor = element_line(color = "black"), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        legend.position = "none",
        text=element_text(size=12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))

# ggsave("./Figures/v4v5_genus.pdf", height = 6, width = 11.4)
# ggsave("./Figures/v3v4_genus.pdf", height = 6, width = 11.4)
ggsave("./Figures/v3v5_genus.pdf", height = 3, width = 10.9)

# Water column
# order_y_by = c("Marine Group II_unclassified","Nitrosopumilaceae_unclassified","Candidatus Nitrosopumilus",
#                "Arctic97B-4 marine group_unclassified", "Rhodopirellula",
#                "Marinimicrobia (SAR406 clade)_unclassified",
#                "Woeseia","SUP05 cluster","SAR86 clade_unclassified",
#                "SAR202 clade_unclassified","SAR324 clade(Marine group B)_unclassified","Polaribacter 1","NS9 marine group_unclassified","NS5 marine group",
#                "Alphaproteobacteria_unclassified","Clade II_unclassified","Clade I_unclassified","Clade Ib","Clade Ia","Subgroup 6_unclassified"),

# Sediment v3v4
# order_y_by = c("Urania-1B-19 marine sediment group","Nitrospira","JG30-KF-CM66_unclassified",
#                "Gemmatimonadaceae_unclassified","BD2-11 terrestrial group_unclassified",
#                "Gammaproteobacteria_unclassified","Woeseia","MBMPE27_unclassified","SAR324 clade(Marine group B)_unclassified","NB1-j_unclassified","bacteriap25_unclassified",
#                "S085_unclassified","SAR202 clade_unclassified",
#                "Mesorhizobium","Alphaproteobacteria_unclassified","Magnetospiraceae_unclassified","Kiloniellaceae_unclassified",
#                "Subgroup 22_unclassified","Subgroup 21_unclassified","Actinomarinales_unclassified"),