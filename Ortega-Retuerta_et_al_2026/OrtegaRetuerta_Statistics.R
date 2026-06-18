################################################################################
# Project: OdiSea Statistical Analysis
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Updated: June 2026
# Description: Script covering basic data processing, fold-change calculations,
#              non-parametric/parametric testing (Wilcoxon, t-test, ANOVA), 
#              and linear correlation modeling.
################################################################################

#####################################
# Load libraries
#####################################
library(readxl)   # For importing modern Excel sheets via read_excel
library(xlsx)     # For alternative Excel workflows via read.xlsx2
library(ggplot2)  # For data visualization and diagnostics plots

##########################
# Load & Filter Data
##########################
# Load main experimental data rows 1-54
Exp <- subset(as.data.frame(read_excel('./OdiSea_data.xlsx', col_names = T)))[c(1:54),]

# Format the 'time' attribute as a categorical factor
Exp$time <- as.factor(Exp$time)

# Order and define explicit levels for the 'treatment' factor grouping
Exp$treatment <- factor(Exp$treatment, levels = c("31.May", "26.Aug"))

# Segment data subsets isolated by seasonal treatments
spring <- subset(Exp, treatment == "31.May")
summer <- subset(Exp, treatment == "26.Aug")

# [Commented Out] Custom subsetting filtering specific treatment and time metrics
# change_spring <- subset(Metadata, 
#                         treatment == "31.May" & time == "t0" | 
#                           treatment == "31.May" & time == "t3")

##########################
# Difference Two Points (Fold-Change & Metrics)
##########################
# Calculate direct fold-change of mean Bacterial Production (BP) at t0 between summer and spring
foldChange = mean(na.omit(t0_summer$BP_nmolCleuL1h1))/mean(na.omit(t0_spring$BP_nmolCleuL1h1))
foldChange

# Calculate percentage decrease of total dissolved amino acids (daa) between spring and summer at t0
percentage = ((mean(na.omit(t0_spring$Total_daa_umolL)) - 
                 mean(na.omit(t0_summer$Total_daa_umolL))) / 
                mean(na.omit(t0_spring$Total_daa_umolL))) * 100
percentage

##########################
# Hypothesis Testing (t-test / Wilcoxon Non-Parametric Framework)
# Wilcoxon test - data are not normally distributed. Sample size should be at least 6!
##########################
# Diagnostic plot checking total amino acids (Total_aa) distribution across treatments at t0
ggplot() + geom_point(data = t0, aes(x=treatment, y = Total_aa,  color = treatment)) 

# Parametric Student's t-test comparing bacterial cells by treatment (Welch's correction enabled via var.equal=FALSE)
t.test(Bacteria_cellsml ~ treatment, data = t0, exact = FALSE, alternative="two.sided", var.equal=FALSE, conf.level=0.95)

# Non-parametric Wilcoxon Rank Sum test assessing bacterial cell discrepancies across treatments
wilcox.test(Bacteria_cellsml ~ treatment, data = t0, exact = FALSE, alternative="two.sided", var.equal=FALSE, conf.level=0.95)

##########################
# Analysis of Variance (ANOVA) & Post-Hoc Controls
##########################
# One-Way ANOVA: Testing Bacterial Production variance against 'day' tracker within Spring group
anova <- aov(BP_ugCLd ~ day, data = spring)
summary(anova)
TukeyHSD(aov(anova)) # Tukey Honest Significant Difference post-hoc pairs test

# One-Way ANOVA: Testing Bacterial Production variance against 'time' factor within Spring group
anova <- aov(BP_ugCLd ~ time, data = spring)
summary(anova)
TukeyHSD(aov(anova))

# Two-Way ANOVA with Interaction: Evaluating main effects and joint effects (time * treatment) on total dataset
anova <- aov(BP_ugCLd ~ time*treatment, data = Exp)
summary(anova)
TukeyHSD(aov(anova))

##########################
# Regression and Correlation
# Reference Model Concept: https://www.geeksforgeeks.org/autocorrelation-and-partial-autocorrelation/?ref=ml_lbp
##########################
# Scatter plot checking Dissolved Inorganic Phosphorus (DIP) over days, overlaid with a linear trend model
ggplot(data = Exp, aes(x =Day , y= DIP))+
  geom_point(aes(color = treatment))+geom_smooth(method = "lm", se = FALSE, color = "black")

# Linear model (regression) evaluating DIP as a function of leucine-incorporated Bacterial Production
lm(DIP ~ BP_pmolCleuL1h1, data = spring)

# [Commented Out] Standard Spearman correlation check between bacteria cell counts and DIP
# cor(Metadata$Bacteria_105cells, Metadata$DIP, method = "spearman", use = "complete.obs")

# Spearman's Rank Correlation test evaluating relationship of DIP against Alkaline Phosphatase Activity (APA) in Spring
cor.test(spring$DIP, spring$APA_uM_h1,  method = "spearman", na.action = "na.exclude")

# Spearman's Rank Correlation test evaluating relationship of DIP against Alkaline Phosphatase Activity (APA) in Summer
cor.test(summer$DIP, summer$APA_uM_h1,  method = "spearman", na.action = "na.exclude")
