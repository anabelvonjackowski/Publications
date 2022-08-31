
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## Correlation 16S and Physiochemical
############################################################################################

library(microbiomeSeq); packageVersion("microbiomeSeq")

# Step: BAC.otu ~ TEP.Area + Galactose + Fucose + CSP.Area + Rhamnose +      Tyrosine + Galactosamine + Glycine
pseq.abs_subset <- subset_samples(pseq.abs, Layer == "100m")
pseq.abs_subset2 <- subset_samples(pseq.abs_subset, Season == "Fall")

physeq.corr <- taxa_level(pseq.abs_subset2, "Genus")

tax_env_cor <- taxa.env.correlation(physeq.corr, "Season", 
                                    # method = "pearman",
                                    method = "spearman",
                     pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 1,
                     num.taxa = 150
                     # ,select.variables = c("Fucose","Galactose" ,"Glycine","Rhamnose" ,"Tyrosine","TEP.Area","CSP.Area")
                     ,select.variables = c("Fucose","Rhamnose","Arabinose","Galactose","Glucose","Mannose_Xylose",
                                           "Galactosamine","Glucosamine","Galacturonic.Acid","Glucuronic.Acid",
                                           "Glycine","Threonine","Arginine","Tyrosine","Valine","Isoleucine","Phenylalanine","Leucine",
                                           "Aspartic.Acid","Glutamic.Acid","Serine","Alanine","GABA")
                     )

# tax_env_cor$Type <- factor(tax_env_cor$Type, levels=c("Summer","Fall"))

## Arrange Environmental Variables
tax_env_cor$Env <- factor(tax_env_cor$Env, levels=c("Threonine","Alanine","Glycine","Glucose","GABA", "Arabinose","Galactose","Galactosamine","Glucosamine","Rhamnose","Aspartic.Acid","Glucuronic.Acid",
                                                    "Fucose","Galacturonic.Acid","Phenylalanine","Leucine","Serine","Mannose_Xylose","Arginine","Tyrosine","Valine","Isoleucine","Glutamic.Acid"))

# tax_env_cor$Env <- factor(tax_env_cor$Env, levels=c("Arabinose","Fucose","Galactose","Glucose","Mannose_Xylose","Rhamnose","Galacturonic.Acid","Glucuronic.Acid","Galactosamine","Glucosamine",
#                                                     "Arginine","Glycine","Isoleucine","Leucine","Phenylalanine","Threonine","Tyrosine","Valine","Alanine","Aspartic.Acid","GABA","Glutamic.Acid","Serine"))

################ summer photic ###################
# tax_env_cor$Taxa <- factor(tax_env_cor$Taxa, levels=c("SAR86_clade_unclassified","SAR11_Clade_Ia","Luteolibacter","NS9_marine_group_unclassified","SUP05_cluster",
#                                                       "NS4_marine_group","SAR11_Clade_II_unclassified","OCS116_clade_unclassified","Marine_Group_II_unclassified",
#                                                       "Candidatus_Actinomarina","OM43_clade","Polaribacter_1", "Aurantivirga","SAR92_clade","Nitrincolaceae_unclassified",
#                                                       "OM60.NOR5._clade","Marinimicrobia_.SAR406_clade._unclassified","Candidatus_Nitrosopumilus","Lentimonas",
#                                                       "Cryomorphaceae_unclassified","Ulvibacter","Formosa","Amylibacter","Ascidiaceihabitans","NS5_marine_group",
#                                                       "Synechococcus_CC9902","SAR11_Clade_IV_unclassified","Planktomarina","SAR116_clade_unclassified","Flavobacteriaceae_unclassified"))
# tax_env_cor$Env <- factor(tax_env_cor$Env, levels=c("Mannose_Xylose","Tyrosine","Alanine","Glycine","GABA",
#                                                     "Galacturonic.Acid","Arabinose","Galactose","Isoleucine",
#                                                     "Phenylalanine","Leucine","Glucuronic.Acid",
#                                                     "Arginine","Valine","Aspartic.Acid","Serine","Rhamnose",
#                                                     "Glutamic.Acid","Glucose","Fucose","Glucosamine","Galactosamine","Threonine"))
########### summer 100 ###################
# tax_env_cor$Taxa <- factor(tax_env_cor$Taxa, levels=c("Pseudohongiella","SAR11_Clade_Ia","SAR86_clade_unclassified","Luteolibacter","SAR11_Clade_Ib","SUP05_cluster","NS9_marine_group_unclassified",
#                                                       "SAR11_clade_unclassified","SAR11_Clade_II_unclassified","OM43_clade","NS4_marine_group","Roseibacillus","OCS116_clade_unclassified",
#                                                       "AEGEAN.169_marine_group_unclassified","Marine_Group_II_unclassified","Candidatus_Nitrosopelagicus","Marinimicrobia_.SAR406_clade._unclassified",
#                                                       "LS.NOB","SAR324_clade.Marine_group_B._unclassified","Candidatus_Nitrosopumilus","Lentimonas","Cryomorphaceae_unclassified","NS5_marine_group",
#                                                       "SAR11_Clade_IV_unclassified","Planktomarina","Flavobacteriaceae_unclassified","Pseudoalteromonas","Aurantivirga","SAR92_clade","Nitrincolaceae_unclassified"))
# 
# tax_env_cor$Env <- factor(tax_env_cor$Env, levels=c("Galacturonic.Acid","Isoleucine","Phenylalanine","Leucine","Glucuronic.Acid","Valine","Aspartic.Acid","Serine","Galactose","Rhamnose",
                                                    # "Arabinose","Arginine","Glutamic.Acid","Mannose_Xylose","Tyrosine","Alanine","Glycine","GABA","Glucose","Threonine", "Fucose","Galactosamine",
                                                    # "Glucosamine"))

########### fall photic ###################
# tax_env_cor$Taxa <- factor(tax_env_cor$Taxa, levels=c("SAR116_clade_unclassified","Amylibacter","Lentimonas","Polaribacter_1","Aurantivirga","SAR11_Clade_IV_unclassified",
#                                                       "Candidatus_Actinomarina","SAR11_Clade_Ia","NS4_marine_group","NS9_marine_group_unclassified","OCS116_clade_unclassified",
#                                                       "OM43_clade", "Nitrincolaceae_unclassified","SAR92_clade","SAR86_clade_unclassified","SAR11_Clade_II_unclassified",
#                                                       "Formosa","Cryomorphaceae_unclassified","Synechococcus_CC9902","Ascidiaceihabitans","NS5_marine_group",
#                                                       "Planktomarina","Ulvibacter","Flavobacteriaceae_unclassified","OM60.NOR5._clade","Marinimicrobia_.SAR406_clade._unclassified",
#                                                       "Marine_Group_II_unclassified","SUP05_cluster","Candidatus_Nitrosopumilus","Luteolibacter"))
# tax_env_cor$Env <- factor(tax_env_cor$Env, levels=c("Threonine","Alanine","GABA","Glucose","Glycine","Arabinose","Galactose","Glucosamine","Rhamnose","Aspartic.Acid",
#                                                     "Glucuronic.Acid","Fucose","Galacturonic.Acid","Phenylalanine","Galactosamine","Leucine","Serine","Isoleucine","Glutamic.Acid","Mannose_Xylose",
#                                                     "Valine","Arginine","Tyrosine"))
########### fall 100 ###################
# tax_env_cor$Taxa <- factor(tax_env_cor$Taxa, levels=c("SAR11_clade_unclassified","OM43_clade","OCS116_clade_unclassified","Nitrincolaceae_unclassified","SAR92_clade","SAR86_clade_unclassified",
#                                                       "SAR11_Clade_II_unclassified","Aurantivirga","Pseudoalteromonas","SAR324_clade.Marine_group_B._unclassified","Candidatus_Nitrosopelagicus",
#                                                       "SAR11_Clade_Ib","SUP05_cluster","LS.NOB","Candidatus_Nitrosopumilus","Luteolibacter","Pseudohongiella","Marine_Group_II_unclassified",
#                                                       "AEGEAN.169_marine_group_unclassified","Marinimicrobia_.SAR406_clade._unclassified","Planktomarina","NS5_marine_group","Cryomorphaceae_unclassified",
#                                                       "Lentimonas","SAR11_Clade_IV_unclassified","SAR11_Clade_Ia","NS4_marine_group","NS9_marine_group_unclassified","Roseibacillus","Flavobacteriaceae_unclassified"))

## Original Formula
# plot_taxa_env(tax_env_cor, aes(x = Env, y = Taxa, fill = Correlation))

## Adapted Forumula
# plot_taxa_env2 <- function(df,...){
#   p <-ggplot2::ggplot(aes(x=Env, y=Taxa, fill=Correlation), data=df)
#   p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")
#   p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
#   p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
#   p<-p+ ggplot2::facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")
#   p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
#   return(p)
# }
# suppressPackageStartupMessages({library(ggplot2)})

plot_taxa_env2(aes(x = Env, y = Taxa, fill = Correlation), df=tax_env_cor)+
  # xlab("Environmental Variables")+
  scale_y_discrete(position = "right")+
  scale_fill_gradientn(colours = c("deepskyblue4","#67a9cf","whitesmoke", "#ef8a62","#B02363"),limits=c(-1,1))+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black"), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))

# write.csv(tax_env_cor, "./Tax_env_summer_photic.csv")
ggsave("./Figures/Correlation_Abundance_Env_fall_100m.pdf")
