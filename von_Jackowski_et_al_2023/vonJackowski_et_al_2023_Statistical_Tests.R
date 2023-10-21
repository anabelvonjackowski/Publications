# © This script was put together by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
setwd("~/.../")
#####################################
#Load libraries
#####################################
library(graphics)
library(phyloseq)
library(vegan); packageVersion("vegan")
library(BiocGenerics)
library(stats)
load(".../x.RData")
#####################################
#1 - OMNIBUS Test for Environmental Parameters #https://www.r-bloggers.com/r-tutorial-series-two-way-omnibus-anova/
#####################################
sediment <- read.csv("./PS101_station_Info_Sediment.csv", header = TRUE, check.names=FALSE)
dataTwoWay <- subset(sediment, select = c(Location, Mean_Depth_cm, Depth_bsf_cm, Chlorophyll, Phaeopigment,  CPE, Contribution,   TOC , TON ,  C_N))
#dataTwoWay$Mean_Depth_cm <- as.character(as.numeric(dataTwoWay$Mean_Depth_cm)) #number --> Charater
dataTwoWay$Mean_Depth_cm <- as.factor(as.numeric(dataTwoWay$Mean_Depth_cm)) #number --> Factor
dataTwoWay <- na.omit(dataTwoWay)
str(dataTwoWay)
anova(lm(Chlorophyll~ Location*Mean_Depth_cm, dataTwoWay)) #two-way Omnibus
anova(lm(Chlorophyll~ Location, dataTwoWay)) #one-way Omnibus

# # eg
# Response: Chlorophyll
#                         Df   Sum Sq   Mean Sq F value   Pr(>F)    
# Mean_Depth_cm           6 0.116456 0.0194093  8.7511 5.14e-07 ***
#   Location              7 0.023202 0.0033145  1.4944   0.1847    
# Mean_Depth_cm:Location 27 0.042664 0.0015801  0.7124   0.8345    
# Residuals              66 0.146383 0.0022179                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# The output of our ANOVA test indicates that the difference between our Mean_Depth means 
# is statistically significant (p < .000) and that the difference between Location is not (p = .18). 
#####################################
##ANOVAs - CPE-chlorophyll-contribution
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r#check-anova-assumptions-test-validity
#####################################
# freshness <- subset(sediment, select = c(Location, Mean_Depth_cm, CPE, Chlorophyll, Contribution))
freshness <- sediment
freshness <- na.omit(freshness)
freshness_0.5 <- freshness[grep("0.5", freshness$Mean_Depth_cm), ]
freshness_1.5 <- freshness[grep("1.5", freshness$Mean_Depth_cm), ]
freshness_2.5 <- freshness[grep("2.5", freshness$Mean_Depth_cm), ]
freshness_4.0 <- freshness[grep("4", freshness$Mean_Depth_cm), ]
freshness_15.0 <- freshness[grep("15", freshness$Mean_Depth_cm), ]
res.CPE <- aov(Contribution ~ Location*Mean_Depth_cm , data = freshness)
res.CPE <- aov(Contribution ~ Location , data = freshness_0.5)
summary(res.CPE)
plot(res.CPE,1) #Homogeneity of variances=numbers indicate outliers

#means
tapply(freshness$Chlorophyll, freshness$Mean_Depth_cm, mean) #yey
#The variances
tapply(freshness$Chlorophyll, freshness$Mean_Depth_cm, var)
#sample size
tapply(freshness$Chlorophyll, freshness$Mean_Depth_cm, length)
#boxplot
boxplot(freshness$Chlorophyll ~ freshness$Mean_Depth_cm)


res.CPE <- aov(C_N ~ Location*Mean_Depth_cm , data = freshness)
res.CPE <- aov(C_N ~ Location , data = freshness_15.0)
summary(res.CPE)

##Tukey
TukeyHSD(res.CPE)

leveneTest(Contribution ~ Location , data = freshness_0.5)
#From the output above we can see that the p-value is not less than the 
#significance level of 0.05. This means that there is no evidence to suggest 
#that the variance across groups is statistically significantly different. 
#Therefore, we can assume the homogeneity of variances in the different 
#treatment groups.

#if leveneTest is violated do:
oneway.test(CPE ~ Location , data = freshness_0.5)
pairwise.t.test(TOC_TON_0.5$C_N, TOC_TON_0.5$Location,
                p.adjust.method = "BH", pool.sd = FALSE)

#check for normality 
plot(res.CPE,2) #numbers indicate normality if they follow the line
#The conclusion above, is supported by the Shapiro-Wilk test on the ANOVA 
#residuals which finds no indication that normality is violated.
aov_residuals <- residuals(object = res.aov2 )
shapiro.test(x = aov_residuals )

#ANOVA assumptions not met? --> kruskal.test
kruskal.test(weight ~ group, data = my_data)

#####################################
#Tukey test or T-Test --> on ANOVA
#####################################
res.CPE <- aov(C_N ~ Location , data = sediment)
TukeyHSD(res.CPE)

or ###t.test
pairwise.t.test(sediment$C_N, sediment$Location,p.adjust.method = "BH") #Benjamini-Hochberg method

# Polaris Vent Reference Seamount Peak 3 Seamount Peak 2 Seamount Saddle Seamount Peak 1N Seamount Peak 1S
# Reference        0.85         -         -               -               -               -                -               
#   Seamount Peak 3  0.90         0.85      -               -               -               -                -               
#   Seamount Peak 2  0.87         0.85      0.90            -               -               -                -               
#   Seamount Saddle  0.79         0.87      0.79            0.79            -               -                -               
#   Seamount Peak 1N 0.79         0.87      0.79            0.79            0.97            -                -               
#   Seamount Peak 1S 0.79         0.85      0.79            0.79            0.87            0.87             -               
#   Southern Slope   0.87         0.87      0.87            0.85            0.85            0.85             0.79            

## Greater than .05 is not statistially significantly difference
## Less or equalt to .05 does have statistically significant difference

#####################################
#ANOVA for rarefaction curves
#####################################
#Omnibus Anova
alpha.diversity <- estimate_richness(all_data, measures = c("Observed","Chao1", "Shannon"))
head(alpha.diversity) # overview
# Observed    Chao1  se.chao1  Shannon
# X10     2746 3790.600  90.03501 5.836359
# X11     2248 2846.962  56.93891 4.097020
# X12     2158 3078.192  92.42434 5.302194
# X13     2886 4139.252 107.57717 5.640149
# X14     2265 3062.942  80.35136 5.633490
# X15     2883 3776.013  78.87574 5.633001

# calculate Anova
data <- cbind(sample_data(all_data), alpha.diversity)
alpha.anova <- aov(Observed ~ Layer*Location, data)
alpha.anova <- aov(Observed ~ Depth, data)
summary(alpha.anova)
# Df Sum Sq Mean Sq F value Pr(>F)  
# Depth_s      2 0.1321 0.06604   4.066 0.0425 *
#   Residuals   13 0.2112 0.01624                 
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Now on subset if previous is significant
# subset your data
all.data.sub<- subset_samples(all_data, Depth == "Epipelagic") #for water
all.data.sub<- subset_samples(all_data, Depth != "Epipelagic") #for water
all.data.sub<- subset_samples(all_data, Layer == "0-1cm") #for sediment
all.data.sub<- subset_samples(all_data, Location != "Hausgarten")

# calculate Anova
alpha.diversity_2 <- estimate_richness(all.data.sub, measures = c("Observed","Chao1", "Shannon"))
head(alpha.diversity_2) # overview
data_2 <- cbind(sample_data(all.data.sub), alpha.diversity_2)

alpha.anova_2 <- aov(Observed ~ General, data_2)
summary(alpha.anova_2)


#####################################
###ANOSIM - Ranking of dissimilarity in ANOSIM and NMDS (non-metric multidimensional scaling) go hand in hand
all.data.sub<- subset_samples(all_data.clr)
all.data.sub<- subset_samples(all.data.sub, Depth_s != "Bottom Water") #for water
all.data.sub<- subset_samples(all.data.sub, Location_s != "Polaris Vent") #for water
all.data.sub<- subset_samples(all_data_MUC, Layer == "Surface (0-1cm)") #for sediment
ano_groups <- get_variable(all.data.sub, "Layer")
ano_groups <- get_variable(all.data.sub, "Depth")
ano_groups <- get_variable(all.data.sub, "Location")
ano_groups <- get_variable(all.data.sub, "Substrate")

unw_pcoa_dm <- phyloseq::distance(all.data.sub, "euclidean")

X1_ano <- anosim(unw_pcoa_dm, ano_groups)
X1_ano$signif
X1_ano$statistic

#####################################
#ANOVA using Adonis function
#####################################
# adonis(data_bray ~ Depth*Location, data = sampledf)

# Calculate bray curtis distance matrix
data_bray <- phyloseq::distance(all_data.clr, method = "euclidean")

# make a data frame from the sample_data
sampledf <- data.frame(all.data.sub@sam_data)

# Adonis test
adonis2(data_bray ~ Depth, data = sampledf)
res.CPE <- aov(data_bray ~ Depth , data = sampledf)
TukeyHSD(res.CPE)
# Homogeneity of dispersion test
beta <- betadisper(data_bray, sampledf$Location)
permutest(beta)
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups     6 0.34161 0.056935 3.8103    999  0.003 **
# Residuals 69 1.03103 0.014942    

#####################################
#PERMANOVA
#hypothesis: Layers we collected samples from have different centroids
set.seed(1234)
#####################################
# Calculate bray curtis distance matrix
# Adonis test
# all.data.sub<- subset_samples(all_data_MUC)
# all.data.sub<- subset_samples(all_data, Location != "Hausgarten") #for water
sampledf <- data.frame(all_data.clr@sam_data)
data_bray <- phyloseq::distance(all_data.clr, method = "euclidean")
adonis2(data_bray ~ Depth*Location, data = sampledf,  permutations=999) #nested = strata=sampledf$Location,
adonis2(data_bray ~ Layer%in%Substrate, data = sampledf,  permutations=999) #nested = strata=sampledf$Location,

# Homogeneity of dispersion test
beta <- betadisper(data_bray, sampledf$Substrate)
# beta <- betadisper(data_bray, sampledf$Depth_s)
permutest(beta)
# Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups     2 0.09918 0.049591 5.9244    999  0.005 **
#   Residuals 73 0.61106 0.008371
#ANSWER is: Yes they have different centroids

#########################
#Chi-2 Test
#####################################
all.data.sub<- subset_samples(all_data, Depth_s == "Epipelagic") #for water
data_bray <- phyloseq::distance(all.data.sub, method = "bray")
sampledf <- data.frame(all.data.sub@sam_data[,1:8])
res.aov <- aov(data_bray ~ Location_r, data = sampledf)
chisq.test(data_bray) 