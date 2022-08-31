
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## LogFC Analysis
############################################################################################

library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2") #BiocManager::install("DESeq2")
library("phyloseq"); packageVersion("phyloseq") #BiocManager::install("phyloseq")
library("cowplot"); packageVersion("cowplot")
library("gage"); packageVersion("gage") #BiocManager::install("gage")
library("dplyr"); packageVersion("dplyr")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')

############################################################################################
###  Calculate Differences in bacterial community composition
############################################################################################
BAC_pruned <- prune_taxa(taxa_sums(pseq.abs)>0,pseq.abs)
BAC <- prune_taxa(taxa_sums(BAC_pruned)>0,BAC_pruned)

# define alpha
alpha_val <- 0.05
q_val <- 0.05

# run DEseq2

BAC.ddsMat <- phyloseq_to_deseq2(BAC, ~Season)
varianceStabilizingTransformation(BAC.ddsMat, blind = TRUE, fitType = "parametric")
BAC.ddsMat <- estimateSizeFactors(BAC.ddsMat)
BAC.ddsMat <- estimateDispersions(BAC.ddsMat)
#geoMeans <- apply(counts(ddsMat), 1, gm_mean)
BAC.DEseq <- DESeq(BAC.ddsMat, fitType="parametric")
BAC.DEseq.res <- results(BAC.DEseq, alpha = 0.05)

############################################################################################
###  Extract Log2FC taxonomic level
############################################################################################
BAC_tax <- as.data.frame(tax_table(BAC_pruned))
BAC_tax$OTU <- rownames(tax_table(BAC_pruned))

BAC_tax %>%
  mutate_if(is.factor, as.character) -> BAC_tax

# Extract the desired taxonomic level
Genera <- unique(BAC_tax$Genus)
colnames(BAC_tax)[6] <- "id"

# Generate list of OTU for each Taxa
OTU.gs <- list()

for (s in 1:length(Genera)){
  n <- Genera[s]
  Order_OTU <- subset(BAC_tax, id == n)
  OTU.gs[[n]] <- Order_OTU$OTU
}


deseq2.fc <- BAC.DEseq.res$log2FoldChange
names(deseq2.fc) <- rownames(BAC.DEseq.res)
exp.fc=deseq2.fc

fc.Order.p <- gage(exp.fc, gsets = OTU.gs, same.dir=TRUE, ref = NULL, samp = NULL)

############################################################################################
###  Plot
############################################################################################

enrch <- rbind(fc.Order.p$greater,fc.Order.p$less)
enrch <-  data.frame(enrch[enrch[,"q.val"]<0.05 &
                                   !is.na(enrch[,"q.val"]),])
enrch$id <- rownames(enrch)
enrch <- unique(merge(enrch,BAC_tax[,c("id","Class")]))

ggplot(data=enrch, aes(y=stat.mean , x=id, label = set.size))+ 
  geom_text(size = 4, aes(y=stat.mean , x=id), nudge_y= -1.4, nudge_x= 0)+
  ylab("Mean log2foldchange")+ 
  geom_point(size = 5, aes(colour = Class))+
  # scale_colour_manual(values = tol21rainbow)+ 
  ylim(-10,10)+
  scale_x_discrete("Family",expand = waiver())+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme(legend.position = "bottom")+
  coord_flip()
