############################################################################################
# MOLA Correlation Environmental & Microbial Co-occurrence Analysis Pipeline
# Author: Anabel von Jackowski
# Framework: microeco & file2meco R6 Class downstream microbiome workflows
# Reference: https://chiliubio.github.io/microeco/
# Author: Anabel von Jackowski (anabelvonjackowski@gmail.com)
# Last Update: June 2024
############################################################################################


# =====================================================
# Load required packages
# =====================================================
library(microeco); packageVersion("microeco")        # Microbial ecology analysis framework
library(file2meco)                                   # Import phyloseq objects into microeco
library(pheatmap)                                    # Heatmap visualization
library(mice)                                        # Missing data handling


# =====================================================
# Prepare 5 m free-living (FL) surface community dataset
# =====================================================

# Subset all samples from 5 m depth
all_data_surface <- phyloseq::subset_samples(
  pseq.abs,
  Depth_m == "5 m"
)

# Select free-living fraction (FL) at 5 m depth
all_data_subset <- phyloseq::subset_samples(
  pseq.abs,
  Depth_m == "5 m" & Filter.size == "FL"
)

# Convert phyloseq object to microeco format
micro <- phyloseq2meco(all_data_subset)

# Extract abundance, taxonomy and metadata tables
otu_table_16S <- micro$otu_table
taxonomy_table_16S <- micro$tax_table
taxonomy_table_16S_tidy <- taxonomy_table_16S %>% tidy_taxonomy

# Separate sample metadata and environmental variables
sample_info_16S <- micro$sample_table[, c(1:11)]
env_data_16S <- micro$sample_table[, -c(1:11)]

# Create microtable object
dataset <- microtable$new(
  sample_table = sample_info_16S,
  otu_table = otu_table_16S,
  tax_table = taxonomy_table_16S
)

# Remove chloroplast and mitochondrial sequences
dataset$filter_pollution(
  taxa = c("mitochondria", "chloroplast")
)

# Rarefy all samples to an equal sequencing depth
dataset$rarefy_samples(sample.size = 10000)

# =====================================================
# Pearson correlation analysis between bacterial genera
# and environmental variables
# =====================================================
#
# Correlations are calculated using Pearson's method.
# P-values are adjusted using the false discovery rate
# (FDR) procedure to account for multiple comparisons.
#
# Correlation matrix:
#     rows    = bacterial genera
#     columns = environmental variables
#
# Significant associations are visualized as a heatmap.
# =====================================================

t4 <- trans_env$new(
  dataset = dataset,
  add_data = env_data_16S[, c(
    "Temperature_C", "Salinity_psu",
    "Chlorophytes", "Cryptophytes",
    "Diatoms", "Dinoflagellates",
    "Haptophytes", "Prasinophytes",
    "Pelagophytes", "Synechococcus",
    "Photoheterotrophic_bacteria_ugL",
    "LAla_nmolL", "LArg_nmolL",
    "LAsp_nmolL", "LGlu_nmolL",
    "LHis_nmolL", "LIleu_nmolL",
    "LLeu_nmolL", "LLys_nmolL",
    "LPhe_nmolL", "LSer_nmolL",
    "LThr_nmolL", "LTyr_nmolL",
    "LVal_nmolL",
    "BAla_nmolL", "DAla_nmolL",
    "DAsp_nmolL", "DGlu_nmolL",
    "DLeu_nmolL", "DSer_nmolL",
    "DVal_nmolL"
  )]
)

# Calculate Pearson correlations with FDR correction
t4$cal_cor(
  use_data = "Genus",
  cor_method = "pearson",
  p_adjust_method = "fdr",
  p_adjust_type = "Env"
)

# Store correlation results
t4$res_cor

# Visualize the correlation matrix as a clustered heatmap
t4$plot_cor(
  pheatmap = TRUE,
  filter_feature = "",
  keep_full_name = TRUE,
  cluster_ggplot = "both",
  color_palette = c(
    colours = c(
      "#2b8cbe",
      "#a6bddb",
      "white",
      "#d8e4bc",
      "#2e8b57"
    )
  )
) +
  theme(
    axis.text.x = element_text(
      size = 12,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor.x = element_blank()
  )

# Export figure
ggsave(
  "./Figures_16S/Correlation_surface_pearson_fdr.pdf",
  height = 6.69291,
  width = 13.7795,
  units = "in",
  dpi = 300
)