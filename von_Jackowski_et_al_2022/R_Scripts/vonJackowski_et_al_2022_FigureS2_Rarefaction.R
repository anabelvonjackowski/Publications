
##############################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## rarefaction and diversity analyses
##############################################################################

# set working directory; save and load
# setwd("/AWI_MPI/FRAM/collaborations/PEBCAO/Rstats")
# load("Pebcao.Rdata")

library(iNEXT) # install.packages("iNEXT")
library(olsrr)
library(cowplot)
library(phyloseq)
library(ggplot2)


############################################################################################
### Calculate Results
############################################################################################
iNEXT <- otu_table(ASV, taxa_are_rows = F)

iNEXT <- iNEXT(as.data.frame(otu_table(iNEXT)), q=c(0),
               datatype="abundance", conf = 0.95, nboot = 100)

## Rarefaction
rarefac <- fortify(iNEXT, type=1)
rarefac.point <- rarefac [which (rarefac$method == "observed"),]
rarefac.line <- rarefac [which (rarefac$method != "observed"),]
rarefac.line$method <- factor (rarefac.line$method,
                               c("interpolated", "extrapolated"),
                               c("interpolation", "extrapolation"))

## Species richness
richness <- ggplot(rarefac, aes(x=x, y=y, colour = site)) +
  geom_line(aes(linetype = method), lwd = 0.5, data = rarefac.line) +
  scale_colour_discrete(guide = FALSE) +
  scale_x_continuous(limits = c(0,1e+5)) +
  labs(x = "Sample size", y = "Species richness") +
  theme_classic(base_size = 12) + theme(legend.position="bottom")

## Sample coverage
cover <- fortify(iNEXT, type=2)
cover.point <- cover [which (cover$method == "observed"),]
cover.line <- cover [which (cover$method != "observed"),]
cover.line$method <- factor (cover.line$method,
                             c("interpolated", "extrapolated"),
                             c("interpolation", "extrapolation"))

ggplot(cover, aes(x=x, y=y, colour = site))+ 
  geom_line(aes(linetype = method), lwd = 0.5, data = cover.line) +
  # scale_colour_discrete(guide = F) +
  scale_colour_manual(values = c("#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24",
                                 "#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24","#df8f24",
                                   "black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black",
                                 "black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black")) +
  scale_x_continuous(limits = c(0,1e+5)) +
  scale_y_continuous(breaks = seq(0.9,1,0.05), limits = c(0.9,1)) +
  labs(x = "Sample size", y = "Sample coverage") +
  theme_classic(base_size = 12) + 
  theme(legend.position="right")

ggsave("./Figures/Rarefaction.pdf",width = 8, height = 8,  bg = "transparent")
