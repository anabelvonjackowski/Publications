
##############################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## basic community analyses
##############################################################################

library(phyloseq)
library(ampvis2); packageVersion("ampvis2")
library(ggplot2)
library(ggConvexHull)
library(dplyr)
library(writexl)
# library(data.table)
source("Amplicon_Colors.R")

# Calculate means
ASV.class<- psmelt(class.rel) %>% as_tibble %>% 
  group_by(Class, Depth, Season, Station) %>%
  summarize_at(c("Abundance"), mean) %>%
  ungroup
ASV.class_1per<- subset(ASV.class, Abundance > 1)

ASV.family <- psmelt(family.rel) %>% as_tibble %>% 
  group_by(Family, Depth, Season) %>%
  summarize_at(c("Abundance"), mean) %>%
  ungroup

ASV.genus <- psmelt(genus.rel) %>% as_tibble %>% 
  group_by(Genus, Depth, Season, Station) %>%
  summarize_at(c("Abundance"), mean) %>%
  ungroup

# Same as plankton
ASV.class$depth <- factor(ASV.class$Depth, levels=c("100m","Below DCM","DCM","Surface"))

ggplot(data = ASV.class_1per,
       aes(x=Station, y=Abundance, fill = Class)) +
  geom_bar(position = "fill", stat = "identity") +
  # ggplot(Data_merged, aes(x=Station, y=Cell.Growth.Rate_d)) +
  facet_grid(Season~Depth)+
  # geom_bar(stat = "identity", position=position_dodge())+
  # scale_x_continuous(trans = "reverse") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(ncol = 4))+
  coord_flip() +
  # ggtitle("Summer (July)")+
  # ggtitle("Fall (Sep/Oct)")+
  ylab("Read Proportions (> 1%)")+
  scale_fill_manual(values = col.class)+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 270, hjust = -.1, vjust = 0.3),
        axis.title=element_text(size=12))

ggsave("./Figures/Bacteria_barplot_class_station_colorblind.pdf", width = 10, height = 8)

############################################################################################
###  Plot Top Taxa using AMPVIS Package
############################################################################################
ampvis$metadata$Depth <- factor(ampvis$metadata$Depth, levels = c("Surface", "DCM","Below DCM","100m"))

ampvis.subset <- amp_subset_samples(ampvis, Layer %in% "Photic")
# ampvis.subset <- amp_subset_samples(ampvis, Layer %in% "100m" & Season %in% "Summer")
ampvis.subset$metadata$Station <- factor(ampvis.subset$metadata$Station, levels = c("N5", "N4","N3","S3","HG4","HG1/2","SV4","SV3","SV2"))

## Order Genera by Class
amp_heatmap(data = ampvis.subset,
  # tax_add = "Class",
 tax_aggregate = "Genus", showRemainingTaxa = TRUE,
 # cluster Photic
 # order_y_by = c("Actinomarinaceae","Cyanobiaceae","Pseudoalteromonadaceae","SAR11_Clade_IV","SAR116_clade","OCS116_clade","Pseudohongiellaceae",
 #                "Thioglobaceae","Saprospiraceae","Cryomorphaceae","Cyclobacteriaceae","Nitrospinaceae","NS9_marine_group","SAR11_Clade_II",
 #                "SAR11_clade_unclassified","AEGEAN-169_marine_group","SAR11_Clade_I","SAR86_clade_unclassified","Nitrosopumilaceae","Puniceicoccaceae",
 #                "Marine_Group_II_unclassified","Marinimicrobia_(SAR406_clade)_unclassified","Nitrincolaceae","Methylophilaceae","Porticoccaceae",
 #                "Rubritaleaceae","Halieaceae","Crocinitomicaceae","Flavobacteriaceae", "Rhodobacteraceae"),
# order_y_by = c("Synechococcus_CC9902","Candidatus_Nitrosopumilus","Marinimicrobia_(SAR406_clade)_unclassified","SAR11_Clade_II_unclassified","Candidatus_Actinomarina",
#                "SUP05_cluster","SAR86_clade_unclassified","OM60(NOR5)_clade","OM43_clade", "NS9_marine_group_unclassified","NS5_marine_group","NS4_marine_group",
#                "SAR116_clade_unclassified","SAR11_Clade_IV_unclassified","SAR11_Clade_Ia","Planktomarina","OCS116_clade_unclassified","Ascidiaceihabitans",
#                "Luteolibacter","Lentimonas","Marine_Group_II_unclassified","SAR92_clade","Nitrincolaceae_unclassified",
#                "Ulvibacter","Polaribacter_1","Formosa", "Flavobacteriaceae_unclassified","Cryomorphaceae_unclassified","Aurantivirga","Amylibacter"),
# #cluster 100m 
# order_y_by = c("Candidatus_Nitrosopumilus","Candidatus_Nitrosopelagicus",
#                "LS-NOB","Marinimicrobia_(SAR406_clade)_unclassified","Pseudoalteromonas",
#                "SAR324_clade(Marine_group_B)_unclassified","SAR11_Clade_II_unclassified","SAR11_Clade_Ib","AEGEAN-169_marine_group_unclassified",
#                "SUP05_cluster","SAR86_clade_unclassified","Pseudohongiella","OM43_clade","NS9_marine_group_unclassified","NS5_marine_group",
#                "NS4_marine_group","Cryomorphaceae_unclassified","SAR11_clade_unclassified","SAR11_Clade_IV_unclassified","SAR11_Clade_Ia",
#                "Planktomarina","OCS116_clade_unclassified","Roseibacillus","Lentimonas", 
#                "Luteolibacter","Marine_Group_II_unclassified","SAR92_clade","Nitrincolaceae_unclassified","Flavobacteriaceae_unclassified","Aurantivirga"),
######  
 tax_show = 30,normalise = T,plot_values = T,
 min_abundance = 0.1,plot_colorscale = "log10",plot_values_size = 4,
 group_by = "Station",facet_by = "Season",
 color_vector = c("deepskyblue4","#67a9cf","whitesmoke", "whitesmoke","#ef8a62","#B02363"))+
  scale_y_discrete(position = "right")+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        legend.position = "none",
        # text=element_text(family="Arial", size=12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggsave("./Figures/Community_heatmap_genus_photic.pdf", height = 6, width = 9)

############################################################################################
### Make textfile
############################################################################################
textmap <- amp_heatmap(ampvis.subset,
                       group_by = "Season",
                       # group_by = c("Season","Station"),
                       tax_aggregate = "Genus",
                       # facet_by = "Season",
                       tax_show = 30, normalise = T,
                       order_y_by = c("SAR11_Clade_Ia","Synechococcus_CC9902","Candidatus_Nitrosopumilus","NS5_marine_group","Marinimicrobia_(SAR406_clade)_unclassified",
                                                     "Candidatus_Actinomarina","SAR11_Clade_II_unclassified","SUP05_cluster","NS9_marine_group_unclassified","SAR86_clade_unclassified",
                                                     "NS4_marine_group","Planktomarina","OCS116_clade_unclassified","SAR116_clade_unclassified",
                                                     "OM60(NOR5)_clade","OM43_clade", "Ascidiaceihabitans", "SAR11_Clade_IV_unclassified","Ulvibacter","SAR92_clade","Polaribacter_1",
                                                     "Marine_Group_II_unclassified","Amylibacter","Luteolibacter","Aurantivirga","Lentimonas",
                                                     "Nitrincolaceae_unclassified","Cryomorphaceae_unclassified","Flavobacteriaceae_unclassified","Formosa"),
                       # order_y_by = c("Candidatus_Nitrosopumilus","SAR11_Clade_Ia","NS9_marine_group_unclassified","SAR86_clade_unclassified","NS5_marine_group",
                       #                               "SAR11_Clade_II_unclassified","Candidatus_Nitrosopelagicus","Marinimicrobia_(SAR406_clade)_unclassified","LS-NOB",
                       #                               "SAR324_clade(Marine_group_B)_unclassified","SAR11_Clade_Ib","AEGEAN-169_marine_group_unclassified","Pseudoalteromonas",
                       #                               "SAR11_Clade_IV_unclassified","Cryomorphaceae_unclassified","Pseudohongiella","SAR11_clade_unclassified","NS4_marine_group",
                       #                               "OM43_clade","Aurantivirga","Roseibacillus","Lentimonas", "Planktomarina","OCS116_clade_unclassified","SAR92_clade",
                       #                               "Flavobacteriaceae_unclassified","Nitrincolaceae_unclassified","Luteolibacter","Marine_Group_II_unclassified","SUP05_cluster"),
                       textmap = TRUE)
textmap
write_xlsx(textmap, './Community_heatmap_genus_photic.xlsx')
