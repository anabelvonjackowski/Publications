################################################################################
# Project: OdiSea Data Visualization
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Updated: June 2026
# Description: Script to load experimental data, preprocess factors, 
#              and generate publication-ready plots using ggplot2.
################################################################################

#####################################
# Load libraries
#####################################
library(readxl)       # For importing modern Excel files (.xlsx) via read_excel
library(xlsx)         # For alternative Excel data manipulation via read.xlsx2
library(ggplot2)      # For data visualization and grammar of graphics
library("wesanderson") # For palette generation inspired by Wes Anderson movies

##########################
# Load & Preprocess Data
##########################
# [Commented Out] Alternative metadata input from "DI_Index" sheet (rows 1-14)
# Metadata <- as.data.frame(read_excel('./aaod21_DI.xlsx', sheet = "DI_Index", col_names = T))[c(1:14),]

# Load main experimental dataset (rows 1-54)
Exp <- subset(as.data.frame(read_excel('./OdiSea_data.xlsx', col_names = T)))[c(1:54),]

# Convert 'time' variable into a categorical factor
Exp$time <- as.factor(Exp$time)

# Convert 'treatment' variable to a factor and strictly order the baseline groups
Exp$treatment <- factor(Exp$treatment, levels = c("31.May", "26.Aug"))

# [Commented Out] Diversity metrics loading (Simpson index subsetting)
# alpha <- read.csv("./richod.csv")
# alpha_subset <- subset(alpha, index == "simpson")

# [Commented Out] Extraction of specific columns for Diversity Index (DI) analysis
# DI <- subset(as.data.frame(read_excel('./OdiSea_data.xlsx', col_names = T, sheet = "DI")))[,c(1:4,56:67)]
# DI$Kaiser_Sum <- as.numeric(DI$Kaiser_Sum)
# head(DI)
# str(DI)

# [Commented Out] Formatting metadata row names and reshaping to long format for analysis
# row.names(Metadata)<- Metadata$ID
# long <- gather(Metadata, key = "Index",value = "DI",
#                # c(12,13,16,18) #absolute
#                c("Kaiser_Sum","PCA_Sum"))

# [Commented Out] Data range verification for the 'Cru' variable
# min(na.omit(Exp$Cru))
# max(na.omit(Exp$Cru))

#####################################
# Data Visualization (ggplot2)
#####################################
ggplot()+
  # [Commented Out] Baseline boxplot for DOC distribution across time/treatment
  # geom_boxplot(Exp, mapping = aes(x=time, y = DOC, color = treatment),  size = 0.5)+
  
  # Scatter plot layer: Standard Dissolved Organic Phosphorus (DOP) data points
  geom_point(Exp, mapping = aes(x=Day, y = DOP, color = treatment),  size = 3, shape = 19)+
  
  # Scatter plot layer: Highlighted outliers using hollow circular borders
  geom_point(Exp, mapping = aes(x=Day, y = DOP_outlier, color = treatment),  size = 3, shape = 21)+
  
  # Trend line layer: Local regression (LOESS) smoothing curve without confidence interval shading
  geom_smooth(Exp, mapping = aes(x=Day, y = DOP, color = treatment), se = FALSE,  span = 0.6)+
  
  # [Commented Out] Alternative axes titles and panel configurations
  # ylab("FDOM (C) Raman Units")+
  # ylab("Alkaline Phosphatase A")+
  # facet_wrap(~Summer, ncol=1)+
  # geom_vline(xintercept=168, color = "blue")+
  # scale_x_discrete(expand = c(0,0))+
  
  # Configure X-axis limits (0 to 35 days) and tick marks at intervals of 5 days
  scale_x_continuous(limits = c(0, 35), breaks = seq(0,35, 5))+
  
  # Configure Y-axis limits
  scale_y_continuous(limits = c(0, 150))+ #, trans="reverse"
  
  # Specify custom shape styles for geometric mapping layers
  scale_shape_manual(values = c(22,23,24))+
  
  # Define precise hex colors for the treatment groups (Brownish-red and Blue)
  scale_color_manual(values = c("#79402E", "#61A8C7FF"))+
  
  # [Commented Out] Statistical summary line and secondary smooth curve overlays
  # stat_summary(fun.y = 'mean', colour = 'grey50', geom = 'line', group = 1)+
  # geom_smooth(se=FALSE, size=1.5, colour= "darkred")+
  
  # Apply clean canvas theme (removes grey gridlines)
  theme_classic()+
  
  # Fine-tune plot styling adjustments for transparent exports
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent", size = 2),
        legend.position = "none", # Hide color legends from the final figure layout
        legend.title = element_text(size=12), 
        # legend.title = element_text(size=25, hjust = 6), 
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent", color = "transparent"), #legend bg
        legend.text = element_text(size=12),
        legend.key.size = unit(1, "cm"),
        axis.text=element_text(size=12),
        # text=element_text(family="Palatino Linotype"),
        # axis.text.x=element_text(angle=0, hjust = 1),
        axis.title=element_text(size=12))

###############################
# Export & Save Figures
###############################
# Export final high-resolution DOP graph
ggsave(filename="./Figures/Odisea_Exp_DOP.pdf", dpi = 300,height = 105, width = 85.3, units = "mm")
