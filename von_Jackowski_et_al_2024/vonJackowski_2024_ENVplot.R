##########################
# Environmental MOLA Data Time-Series Processing and Visualization
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Last Update: June 2024
##########################

##########################
# Load Packages
##########################
library(scales); packageVersion("scales")       # For axis scaling transformations and numeric labeling formatting
library(ggplot2); packageVersion("ggplot2")     # For structuring data graphics layers
library(readxl)                                 # For reading datasets directly out of Excel worksheets
library(lubridate)                              # For handling calendar transformations (Note: detach command commented out inline)
library(dplyr)                                  # For grouping, calculating pipelines, and data aggregation
library(tidyr)                                  # For structural data pivoting (e.g., transforming wide data to long via gather)

##########################
# Calculate Means
##########################
# Aggregate chlorophyll concentration metadata by month to resolve basic descriptive values
df_mean_std <- ENV_subset %>%
  group_by(Month_) %>%
  summarise_at(vars(Chla_ugL), list(mean=mean, sd=sd)) %>% 
  as.data.frame()
colnames(df_mean_std)

# Note: This assigns values from a function-like call utilizing the original 'ENV' frame to a new frame column
df_mean_std$Month <- df_mean_std(ENV$Month)

##########################
# Load and Structure Datasets
##########################
# Load surface dataset at 5 meters depth filtering for free-living (FL) bacteria fractions
Metadata <- subset(as.data.frame(read_excel('./mola_data.xlsx', sheet = "mola")),Depth_m == "5 m" & Filter.size == "FL")
Metadata$Year <- as.factor(Metadata$Year)
# Metadata$Date <- as.factor(Metadata$Date)

# Note on Excel Date Tracking: Chronological values imported from Excel track as days elapsed since 1899-12-30. 
# Multiplying by 86400 converts days to seconds to parse correctly via POSIXct in Central European Time (CET).
Metadata$Date <- as.POSIXct(as.numeric(Metadata$Date)*86400, origin = "1899-12-30", tz="CET")
Metadata$Month <- factor(Metadata$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# head(Metadata)[,c(1:10)]

# Load historical context data sheet covering the top 125 rows of records
MetadataXL <- as.data.frame(read_excel('./mola_data_profile.xlsx', sheet = "mola20220304"))[c(1:125), ]
MetadataXL$Year <- as.factor(MetadataXL$Year)

# Isolate deep basin control checks captured at 500 meters depth
D500m <- subset(as.data.frame(read_excel('./mola_data.xlsx')),Depth_m == "500 m" & Filter.size == "FL")
D500m$Year <- as.factor(D500m$Year)

# Check data boundaries for dissolved hydrolyzable amino acid variants before plotting layers
min(na.omit(D500m$DDHAA_nmolL))
max(na.omit(D500m$DHAA_umolL))
# Metadata_subset <- subset(Metadata, Sample_ID != "sa2160" & Sample_ID != "sa2156")

# [Commented Out] Alternative Environmental Core dataset configuration logic
# ENV <- as.data.frame(read_excel('./molaenv.xlsx'))
# ENV$Year <- as.factor(ENV$Year)
# ENV$Filter.size <- as.factor(ENV$Filter.size)
# ENV$Date <- ymd(ENV$Date)
# ENV$Month <- factor(ENV$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# 
# ENV_subset <- subset(ENV, Filter.size == "0.2")
# ENV_surface <- subset(ENV, Filter.size == "0.2" & Depth_m == "5 m")
# ENV_deep <- subset(ENV, Filter.size == "0.2" & Depth_m == "500 m")
# colnames(ENV_subset)

# Simple Overview Plotting Layer
ggplot()+
  # Map 500m deep data points using Julian dates on the horizontal axis
  geom_point(D500m, mapping = aes(x=JulianDate, y = DDHAA_nmolL, shape = Year), fill = "black", size = 3)+
  
  # Fit local regression smoothing curve mapping over deep target metrics
  geom_smooth(D500m, mapping = aes(x=JulianDate, y = DDHAA_nmolL), size = 0.5,  se = FALSE, color = "black")+
  # geom_point(Metadata, mapping = aes(x=JulianDate, y = LNA_cellsmL, shape = Year), fill = "#3B9BB3", size = 3)+
  # geom_smooth(Metadata, mapping = aes(x=JulianDate, y = LNA_cellsmL), size = 0.5,  se = FALSE, color = "#3B9BB3")+
  # geom_point(Metadata, mapping = aes(x=JulianDate, y = DOC_umolL, shape = Year), fill = "black", size = 3)+
  # geom_smooth(D500m, mapping = aes(x=JulianDate, y = DOP_umolL), size = 0.5,  se = FALSE, color = "black")+
  # geom_point(D500m, mapping = aes(x=JulianDate, y = DOC_umolL, shape = Year), fill = "#FFC000", size = 3)+
  # ylab("Bacteria (x10^5 cells mL-1)")+
  # facet_wrap(~Freshwater, ncol=1)+
  # geom_vline(xintercept=168, color = "blue")+
  # geom_vline(xintercept=210, color = "blue")+
  
  # Format full annual temporal grid axis matching start days of calendar months
  scale_x_continuous(expand = c(0,0), limits = c(0, 365), breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335))+
  # scale_y_continuous(expand = c(0,0),limits = c(0, 10e+05),labels = scales::scientific)+
  
  # Explicitly state visual point styles for tracking categorical differences across years
  scale_shape_manual(values = c(22,23,24))+
  # scale_color_manual(values = c("darkblue", "#FFC000", "red"))+
  # stat_summary(fun.y = 'mean', colour = 'grey50', geom = 'line', group = 1)+
  # geom_smooth(se=FALSE, size=1.5, colour= "darkred")+
  theme_classic()+
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent", color = "transparent", size = 2))

###############################
# Save Output Plots
###############################
# Export basic high-resolution graphic matching target environmental directory parameters
ggsave(filename="./Figures_Environmental/MOLA_Bacteria_groups.pdf", height = 3.9, width = 4.25, dpi = 300) 

###############################
# Reshape Wide Matrices to Long Format
###############################
# Note on Data Subsetting and Reshaping Strategy:
# The following sequential gather calls overwrite the 'long' object reference variable repeatedly. 
# In normal workflows, each block corresponds to generating a specific multi-variable line chart below.

# Pivot 1: Reshaping diverse phytoplankton community counts
long <- gather(Metadata, key = "variable",value = "value",
               # c(12,13,16,18) #absolute
               c("Prasinophytes","Chlorophytes", "Cryptophytes", "Diatoms", 
                 "Dinoflagellates", "Haptophytes","Pelagophytes", "Synechococcus", "Photoheterotrophic_bacteria_mgL"))

# Pivot 2: Reshaping total combined and fundamental standalone amino acid values (filtering out outlier labels)
Metadata_subset <- subset(Metadata, Sample_ID != "S240FL" &  Sample_ID != "S168FL" )
long <- gather(Metadata_subset, key = "variable",value = "value",
               c("LDHAA_nmolL","DDHAA_nmolL", "GABA_nmolL", "Gly_nmolL"))

# Pivot 3: Reshaping second grouping tier of L-form amino acid concentrations
long <- gather(Metadata_subset, key = "variable",value = "value",
               c("LSer_nmolL","LTyr_nmolL",
                 "LHis_nmolL","LLys_nmolL","LVal_nmolL"))

# Pivot 4: Reshaping third grouping tier of L-form amino acid targets
long <- gather(Metadata_subset, key = "variable",value = "value",
               c("LAsp_nmolL","LArg_nmolL", "LGlu_nmolL", "LThr_nmolL", 
                 "LAla_nmolL","LPhe_nmolL", "LIleu_nmolL", "LLeu_nmolL"))

# Pivot 5: Reshaping chiral D-form amino acid variants and unique degradation biomarkers
# Note: Trailing comma within vector selection is kept exactly intact to respect structural choices.
long <- gather(Metadata, key = "variable",value = "value",
               c("BAla_nmolL", "DAla_nmolL",
                 "DAsp_nmolL", "DGlu_nmolL","DSer_nmolL", 
                 "DVal_nmolL", "DLeu_nmolL", ))

# Pivot 6: Reshaping high/low nucleic acid (HNA/LNA) active microbial populations
Metadata_subset <- subset(Metadata, Sample_ID !="S238FLa")
long <- gather(Metadata_subset, key = "variable",value = "value",
               c("LNA_cellsmL", 
                 "HNA_cellsmL"))

# Pivot 7: Reshaping 500m baseline deep water hydrolyzable amino acid profiles
D500m_subset <- subset(D500m,Sample_ID !="D168FL" &  Sample_ID != "D240FL" & JulianDate !="99")
long <- gather(D500m, key = "variable",value = "value",
               c("LDHAA_nmolL","DDHAA_nmolL"))

# Complex Multi-Variable Trend Matrix Visualization
ggplot()+
  # geom_point(Metadata, mapping = aes(x=JulianDate, y = PO4_umolL, shape = Year),fill ="#3B9AB2", size = 3)+
  # geom_point(Metadata, mapping = aes(x=JulianDate, y = DOP_umolL, shape = Year),fill ="black", size = 3)+
  # geom_smooth(MetadataXL, mapping = aes(x=JulianDay, y = POC_mmol_m2,color = "black"), size = 0.5,  se = FALSE)+
  
  # Plot points color-mapped and fill-mapped by individual pivoted components
  geom_point(long, mapping = aes(x=JulianDate, y = value, shape = Year, color = variable, fill = variable), size = 4)+
  
  # Render smooth trend trajectories distinct for each reshaped column vector 
  geom_smooth(long, mapping = aes(x=JulianDate, y = value,color = variable), size = 0.5,  se = FALSE)+
  # geom_vline(xintercept=168, color = "blue")+
  # geom_vline(xintercept=210, color = "blue")+
  # geom_hline(yintercept = 37.9, linetype = "dashed")+
  
  # Format horizontal time bounds to capture complete yearly distribution
  scale_x_continuous(expand = c(0,0), limits = c(0, 365), breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335))+
  
  # Strip scientific float rendering style to print whole integers on vertical tracking scales
  scale_y_continuous(expand = c(0,0),limits = c(0, 1200), labels = function(x) sprintf("%.0f", x))+
  # stat_summary(fun.y = 'mean', colour = 'grey50', geom = 'line', group = 1)+
  
  # Color arrays defined matching custom container variables 'my_color_2'
  scale_fill_manual(values = my_color_2)+
  scale_color_manual(values = my_color_2)+
  # scale_fill_manual(values = c("LAla_nmolL"="#0a4968","LArg_nmolL"="#3B9AB2","LAsp_nmolL"="#0B775E","LGlu_nmolL"="#02401B",
  