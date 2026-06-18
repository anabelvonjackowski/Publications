################################################################################
# Project: OdiSea PCA Variable Factor Map
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Updated: June 2026
# References and Guides for this methodology:
# https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
################################################################################

#####################################
# Load libraries
#####################################
library(FactoMineR); packageVersion("FactoMineR") # For advanced exploratory multivariate data analysis
library(tidyverse)    # For data manipulation, tidying, and pipe workflows
library(devtools)     # For package development utilities and installing external packages
library(ggbiplot)     # For generating ggplot2-based biplots from PCA outputs
library(ggplot2)      # For data visualization and precise layout adjustments
library(plyr)         # For data splitting, applying functions, and combining results
library(scales)       # For internal scale mapping, breaks, and label formatting
library(grid)         # For low-level graphical grid engine layouts
library(readxl)       # For fast, dependency-free importing of Excel worksheets
library(writexl)      # For read_excel utilities and light export features

library(xlsx)         # For legacy Excel (.xlsx) workflows requiring read.xlsx2
#####################################
# Load & Preprocess Data
#####################################
# Load raw PCA target data from the specific worksheet "PCA"
data_AA <- as.data.frame(read_excel('./aaod21_DI.xlsx', sheet = "PCA", col_names = T))
# data_AA <- read.xlsx2('./aaod21_DI.xlsx', sheetName = "PCA", row.names=1,  colClasses=NA,
# as.data.frame = TRUE, header = TRUE)

# Assign the 'ID' column as explicit row labels for identification
row.names(data_AA)<- data_AA$ID

# Format the 'Day' track variable as a discrete categorical factor
data_AA$Day <- as.factor(data_AA$Day)

# data_AA$Date <- ymd(data_AA$Date)
# data_AA$Month <- factor(data_AA$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov"))

# data_AA <- subset(data_AA, Sample != "sa53" & Sample != "sa57" & Depth_m == "5 m")

#####################################
# Run PCA Analysis
#####################################
# Isolate descriptive environmental data vectors (Columns 1-3) for group aesthetics
ir.AA <- data_AA[,c(1:3)] 
str(ir.AA)

# Isolate target Amino Acid concentrations [mol%] for the analysis (excluding GABA)
# Includes both standard L- and D- form Amino Acids
data <-subset(data_AA, select = c("Asx" ,"Glx","Ser" ,"Thr" ,"Arg","Ala","Tyr","Val","Phe",
                                  "Leu","Gly","His","Lys","Met"))
colnames(data)
str(data)

# View(data)

# ir.AA$depth <- factor(ir.AA$depth, levels=c("Surface","above DCM","DCM","below DCM","100m"))
# ir.AA$Timing <- factor(ir.AA$Timing, levels=c("Summer","Autumn"))


# Compute Principal Component Analysis (PCA)
# Highly recommended: center = TRUE and scale = TRUE to normalize different variable scales
ir.pca <- prcomp(data,center = TRUE, scale = TRUE) 

# Print the basic standard deviations and rotation summary matrix
print(ir.pca)

# Extract variable loadings (rotation matrix) and convert to a dataframe for exporting
ir.pca$rotation
y=as.data.frame(ir.pca$rotation) 

# Export Results:
# The print method returns the standard deviation of each of the four PCs,
# and their rotation (or loadings), which are the coefficients of the linear
# combinations of the continuous variables
# write.csv(y,'output_PCA_Abb_Rotation.csv', row.names=T, col.names=T)


# Metadata overview:
# center and scale refers to respective mean and standard deviation of the variables that are used for normalization prior to implementing PCA
# The rotation measure provides the principal component loading. Each column of rotation matrix contains the principal component loading vector.
names(ir.pca)

# Extract individual principal component coordinates (score vectors)
ir.pca$x  

x=as.data.frame(ir.pca$x) 

# Export individual object coordinates for backup or alternative plotting
# write.csv(x, 'output_Scores.csv', row.names=T, col.names=T)


#####################################
# Plot Biplot Graph
#####################################
# Simple Biplot Check
# biplot(ir.pca, scale = 0)
# dev.off()
# plot method
# plot(ir.pca, type = "l")
# summary method
# The first row describe again the standard deviation associated with each PC.
# The second row shows the proportion of the variance in the data explained by each 
# component while the third row describe the cumulative proportion of explained variance. 
# summary(ir.pca)


# Verify available groupings before plotting the sophisticated biplot
colnames(ir.AA)

# Build custom publication-ready PCA Biplot via ggbiplot
ggbiplot(ir.pca, 
         choices = 1:2,       # Plot PC1 against PC2
         varname.size = 5,    # Text size of variable labels
         var.axes = T,        # Include directional vector arrows
         # groups = ir.AA$Season,
         # labels = ir.AA$W_Layer,
         obs.scale = 1, var.scale = 1
) +  
  # scale_color_manual(name="dCCHO [mol%]",values=c("navy","orange","#72FAA2","#A8FA72"))+
  
  # Map sample points color to season status and geometric shape to specific time points
  geom_point(aes(color=ir.AA$Summer, shape = ir.AA$Day), size=4)+
  
  # Tighten coordinate axes layouts without padding gaps
  scale_x_continuous(expand=c(0,0), limits = c(-6,4))+
  scale_y_continuous(expand=c(0,0), limits = c(-4,4))+
  # geom_text_repel(aes(label = ir.AA$Month))+
  # scale_color_manual(values = c("Jan"="#3B9AB2","Feb"="#0B775E", "Mar"="#3db229", "Apr"="#F3DF6C", "May"="#DC863B", "Jun"="#6d4418", 
  #                               "Jul"="#C7B19C", "Aug"="#A42820", "Sep"="#E6A0C4", "Oct"="#F8AFA8", "Nov"="black"))+
  
  # Apply distinct colors for Summer groups
  scale_color_manual(values=c("gold","cadetblue3"))+
  # guides(fill=guide_legend(override.aes=list(shape=21)))+
  
  # Configure descriptive titles for the plot legends
  labs(shape="Timepoint", color="Summer")+
  
  # Define explicitly different glyph marker assignments
  scale_shape_manual(values = c(16,17,18,25,8,15)) +
  ggtitle("dHAA [mol%]")+
  
  # Apply white panel backing framework
  theme_bw() +
  
  # Fine-tune theme controls for transparent graphic export properties
  theme(plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        # panel.background = element_rect(fill = "#faf0e6", color = NA), # bg of the plot
        panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        # panel.border = element_rect(fill = "transparent", size = 1),
        strip.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        axis.text=element_text(size=12),
        # axis.text.x = element_text(angle = 90, hjust = 1, size=16),
        axis.title=element_text(size=12))

#####################################
# Save Plot
#####################################
# Export high-resolution landscape PDF vector file matching dimensions (A4 scale)
ggsave(file = "./Odisea_DHAA_PCA.pdf", dpi = 300,height = 210, width = 297, units = "mm")