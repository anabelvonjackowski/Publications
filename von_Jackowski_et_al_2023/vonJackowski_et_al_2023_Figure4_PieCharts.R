# © This script is written by Anabel von Jackowski @avonjackowski

#####################################
#Load
#####################################
library(phyloseq)
library(ampvis2)
library(dplyr)
load("./data_import_v3v4.RData")
load("./data_import_v4v5.RData")
#####################################
#Plot top taxa
#####################################
pseq.abs = filter_taxa(all_data, function(x) sum(x > 3) > (0.03 * length(x)), TRUE)
pseq.rel = transform_sample_counts(pseq.abs, function(x) x / sum(x) * 100) 
class.rel <- tax_glom(pseq.rel, taxrank = "Class")
all_data_df <- data.frame(OTU = rownames(phyloseq::otu_table(pseq.abs)@.Data),
                      phyloseq::otu_table(pseq.abs)@.Data,
                      phyloseq::tax_table(pseq.abs)@.Data,
                      check.names = FALSE)

# Extract metadata from phyloseq; format and combine
all_data_df.env <- data.frame(phyloseq::sample_data(pseq.abs), check.names=F)
all_data_df.env <- cbind(Sample = rownames(all_data_df.env), all_data_df.env) 
ampvis <- amp_load(all_data_df, all_data_df.env)

#Plot top 10 as heatmap
amp_heatmap(ampvis, tax_show = 10, tax_aggregate = "Class",showRemainingTaxa = TRUE,measure = "mean"
            # ,group_by = c("Location", "Layer")
  )
# Alphaproteobacteria  Gammaproteobacteria Bacteroidia Deltaproteobacteria Dehalococcoidia
# v3v4 specifics: Plankcomycetacia Phycisphaerae  Acidimicrobiia Gemmatimonadetes  Nitrospira  Plankcomycetacia
# v4v5 specifics: Nitrosopheria Marinimicrobia (SAR406 clade)_unclassified Thermoplasmata Verrucomicrobiae

# Calculate means
# ASV.class$Class[ASV.class$Abundance < 5] <- "Other Bacteria"

# v3v4 Classes
# ASV.class<- psmelt(class.rel) %>% as_tibble %>% group_by(Class, Layer,Location) %>% summarize_at(c("Abundance"), mean) %>% ungroup
# ASV.class$Class[ASV.class$Class != "Alphaproteobacteria" & ASV.class$Class != "Gammaproteobacteria" & ASV.class$Class != "Dehalococcoidia" & ASV.class$Class != "Acidimicrobiia" & ASV.class$Class != "Deltaproteobacteria" & 
#                 ASV.class$Class != "Phycisphaerae" & ASV.class$Class != "Bacteroidia" & ASV.class$Class != "Gemmatimonadetes" & ASV.class$Class != "Nitrospira" & ASV.class$Class != "Planctomycetacia"] <- "Other Bacteria"
# ASV.class$Layer <- factor(ASV.class$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))
# ASV.class$Location <- factor(ASV.class$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle","Karasik Seamount", "Southern Slope","S-Reference"))

# v4v5 Classes
ASV.class<- psmelt(class.rel) %>% as_tibble %>% group_by(Class, Depth,Location) %>% summarize_at(c("Abundance"), mean) %>% ungroup
ASV.class$Class[ASV.class$Class != "Alphaproteobacteria" & ASV.class$Class != "Gammaproteobacteria" & ASV.class$Class != "Bacteroidia" & ASV.class$Class != "Deltaproteobacteria" & ASV.class$Class != "Dehalococcoidia" &
                  ASV.class$Class != "Planctomycetacia" & ASV.class$Class != "Nitrosopheria" & ASV.class$Class != "Marinimicrobia (SAR406 clade)_unclassified" & ASV.class$Class != "Thermoplasmata" & ASV.class$Class != "Verrucomicrobiae"] <- "Other Bacteria"
ASV.class$Depth <- factor(ASV.class$Depth, levels=c("Surface","200 m", "400 m", "600 m", "Bottom Water", "Bathypelagic"))
ASV.class$Location <- factor(ASV.class$Location, levels=c("Reference Site","Northern Seamount", "Central Seamount","Karasik Seamount"))
unique(ASV.class$Class)

# ASV.class_surface <- subset(ASV.class, Layer == "Surface (0-1cm)")
# ASV.class_suburface <- subset(ASV.class, Layer == "Subsurface (1-5cm)")
# ASV.class_deep <- subset(ASV.class, Layer == "Deep Layer (14-16cm)")

#Plot Pie Charts
ggplot(data=ASV.class, aes(x = "", y = Abundance, fill = Class)) +
  geom_bar(data=ASV.class,aes(y=Abundance, fill=Class, color=Class),
           width=1, stat = "identity") +
  coord_polar("y", start=0) +
  # facet_grid(Layer~Location)+
  facet_grid(Depth~Location)+
  guides(fill = guide_legend(ncol = 2)) + 
  scale_color_manual(values = my_color_10) +
  scale_fill_manual(values = my_color_10) +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position= "right",
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1),
        strip.text.y = element_text(size = 12)
  ) 

ggsave("./Figures/v4v5_Class_Pie.pdf", width = 15, height = 10, dpi = 400)

# save.image(file='./betadiversity_barplots_v3v4.Rdata')
