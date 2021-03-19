# Linear Modell Mixed Model

# The statistical software R (2019) was used to evaluate the data. The data evaluation startet with 
# the definition of an appropriate statistical mixed model (Laird and Ware, 1982; Verbeke and 
# Molenberghs, 2000).  
# The model included Season (Level 1, Level 2) and Tiefe (...), as well as their interaction term 
# as fixed factors. The Station was regarded as random factor. 
# The residuals were assumed to be normally distributed and to be homo/heteroscedastic with 
# respect to the different levels of Season and Tiefe. These assumptions are based on a graphical 
# residual analysis. 
# Based on this model, a Pseudo R^2 was calculated (Nakagawa and Schielzeth, 2013) and an analysis of 
# variances (ANOVA) was conducted, followed by multiple contrast tests (e.g., see 
# Bretz et al., 2011) in order to 
# compare the several levels of the influence factors, respectively. If Season and depth had no 
# significant interaction, then the levels of the respective factor was pooled.

###############################
# Load Libraries
###############################
library(gplots)
library(nlme)
library(piecewiseSEM)
library(multcomp)
library(lsmeans) 
library(car)
library(dplyr)

###############################
# Load Data
###############################
tabelle <- data.merge %>%
  dplyr::select(Season,W_Type,Station,W_Layer,
                Temperature_C, Salinity_PSU,
                Chla_ug_L, DOC_umol_L,Labile.DOC,Labile_DOC_umolL, 
                dCCHO_nmolL,dHAA_nmolL,
                TEP_Area_cm2_L,CSP_Area_cm2_L,
                Bacteria_cells_ml_105,HNA_ml,LNA_ml,BP1_ugCL.1d.1)
str(tabelle)


###############################
# Subset indididual Parameters
###############################
tabneu <- droplevels(subset(tabelle, !is.na(Bacteria_cells_ml_105)))

# windows(10,7); par(mfrow=c(1,1), mar=c(8, 4, 4, 2) + 0.1)
boxplot2(Bacteria_cells_ml_105 ~ Season+W_Type, data=tabneu, main="DOC_micromol_L", xlab="", ylab="", las=3,
  col=rainbow(length(levels(tabneu$Season))))

###############################
# Model and Anova
###############################
# W_Layer or W_Type
mod1 <- lme(Bacteria_cells_ml_105 ~ Season*W_Type,
            data=tabneu, random=~1|Station,
             # weights=varIdent(form=~1|Season*W_Type),
            na.action=na.exclude,
            control=list(maxIter=300,msMaxIter=300,niterEM=300,msMaxEval=300,opt="nlminb"))
mod2 <- lme(Bacteria_cells_ml_105 ~ 0 + Season:W_Type,
            data=tabneu, random=~1|Station,
             # weights=varIdent(form=~1|Season*W_Type),
            na.action=na.exclude,
            control=list(maxIter=300,msMaxIter=300,niterEM=300,msMaxEval=300,opt="nlminb"))
AIC(mod1,mod2)

###############################
# Residual Analysis
###############################
plot(mod1)

summary(mod1)


# Pseudo-R^2
rsquared(mod1)


# ANOVA
anova(mod1)


###############################
# 1.: W_Layer?
###############################

# no significant interaction:
vergl1a <- lsmeans(mod1, specs="W_Type", contr="pairwise")$contrasts
summary(as.glht(vergl1a, by=NULL, alternative="two.sided"))
# 
# significant interaction:
vergl1b <- lsmeans(mod1, specs="W_Type", by=c("Season"), contr="pairwise")$contrasts
summary(as.glht(vergl1b, by=NULL, alternative="two.sided"))

###############################
# 2.: Season?
###############################
# no significant interaction:
vergl1a <- lsmeans(mod1, specs="Season", contr="pairwise")$contrasts
summary(as.glht(vergl1a, by=NULL, alternative="two.sided"))

# significant interaction:
vergl1b <- lsmeans(mod1, specs="Season", by=c("W_Type"), contr="pairwise")$contrasts
summary(as.glht(vergl1b, by=NULL, alternative="two.sided"))

###############################
# Correlation Analysis

# https://stackoverflow.com/questions/19012529/correlation-corrplot-configuration

# 5e-4, "***", = 0.005
# < 5e-3, "**", =< 0.005
# normals:
# < 1e-3, "**", =0.001
# < 1e-2, "*", = 0.01
# < 5e-2, "*", = 0.05
###############################

library(corrgram)
library(ellipse)

tabelle <- data.merge %>%
  dplyr::select(Season,W_Type,
                Temperature_C, Salinity_PSU,
                Chla_ug_L, DOC_umol_L,
                # Labile_DOC_umolL,
                Labile.DOC,
                dCCHO_umolL,#dCCHO_umolCL,
                dHAA_umolL,#dHAA_umolCL,dHAA_umolNL,
                TEP_Area_cm2_L,CSP_Area_cm2_L,
                Bacteria_cells_ml_105,BP1_ugCL.1d.1)

tabelle$Ges <- as.factor(paste(tabelle$W_Layer, tabelle$Season))


str(tabelle)

tabneu <- na.omit(tabelle[,c(3:13,length(colnames(tabelle)))])

tabneu

###############################
# naive Correlation
###############################

nk <- cor(tabneu[,-c(length(colnames(tabneu)))])
round(nk, 3)
# write.csv(nk, file = "./Season_Naive_Correlations.csv")

# pairs(tabneu[,-c(length(colnames(tabneu)))], panel=panel.smooth)
corrgram(tabneu[,-c(length(colnames(tabneu)))],
         order = F,
         lower.panel=NULL,
         upper.panel=panel.shadeNtext,
         col.regions = colorRampPalette(c("navy","#007FFF","#3B9AB2", "#EEEEEE","salmon", "red"))
         # col.regions = colorRampPalette(c("blueviolet","blue3","#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red")) #jetcolors
         )
  

# corrgram(tabneu[,-c(length(colnames(tabneu)))], lower.panel=panel.pts, upper.panel=panel.conf,
  # diag.panel=panel.density)
# windows(10,7); 
corrgram(tabneu[,-c(length(colnames(tabneu)))], upper.panel=panel.pts)

###############################
# Mixed Correlation by Season
###############################

Corrs <- list()
for (i in 1:length(levels(tabneu$Ges))) {
  Corrs[[i]] <- round(cor(droplevels(subset(tabneu,Ges==levels(tabneu$Ges)[i]))[,-c(length(colnames(tabneu)))]),3)
}
names(Corrs) <- levels(tabneu$Ges)
Corrs
# write.csv(Corrs, file = "./Season_Specific_Correlations.csv")

###############################
# Average Mixed Correlations
###############################

mixed <- matrix(0,nrow=nrow(Corrs[[1]]),ncol=ncol(Corrs[[1]]))
for (i in 1:length(levels(tabneu$Ges))) {
  mixed <- mixed+Corrs[[i]]
}
dsk <- mixed/length(levels(tabneu$Ges))
dsk

write.csv(dsk, file = "./Season_Mean_Correlation_WType.csv")
# Differenz der "wei?en" und der Durchschnittekorrelationen:

nk-dsk

# Signifikanzen:
#                 T = r*sqrt(n-2) / sqrt(1-r^2)
#               T^2 = r^2*(n-2) / (1-r^2)
#       T^2*(1-r^2) = r^2*(n-2)
#     T^2 - T^2*r^2 = r^2*(n-2)
#               T^2 = r^2*(n-2) + T^2*r^2
#               T^2 = r^2*(n-2+T^2)
#   T^2 / (n-2+T^2) = r^2
#                 r = T / sqrt((n-2+T^2))
ts <- qt(0.975, df=nrow(na.omit(tabneu))-2, lower.tail=TRUE)
# quantile function = vector of probability, dataframe = 138 rows
rmin <- ts / sqrt(ts^2+nrow(na.omit(tabneu))-2)
# quantile function / sqrt
round(rmin,3)
# -> alle Korrelationen mit |r| > rmin sind signifikant
#    alle nicht-signifikanten hier auf 0 gesetzt

panel.shadeNtext <- function (x, y, corr = NULL, col.regions, ...) 
{
  corr <- cor(x, y, use = "pair")
  results <- cor.test(x, y, alternative = "two.sided")
  est <- results$p.value
  stars <- ifelse(est < 1e-3, "***", 
                  ifelse(est < 1e-2, "**", 
                         ifelse(est < 5e-2, "*", "")))
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                               length = ncol + 1), include.lowest = TRUE))
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
       border = NA)
  box(col = "lightgray")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- formatC(corr, digits = 2, format = "f")
  cex.cor <- .8/strwidth("-X.xx")
  fonts <- ifelse(stars != "", 2,1)
  # option 1: stars:
  text(0.5, 0.4, paste0(r,"\n", stars), cex = cex.cor)
  # option 2: bolding:
  #text(0.5, 0.5, r, cex = cex.cor, font=fonts)
}
