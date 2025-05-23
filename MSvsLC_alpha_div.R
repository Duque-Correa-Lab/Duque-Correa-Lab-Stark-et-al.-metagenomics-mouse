# Load required libraries
library(vegan)         # For diversity indices (e.g., Shannon)
library(ggplot2)       # For plotting
library(stringr)       # String manipulation
library(dunn.test)     # For post-hoc pairwise testing
library(dplyr)         # Data manipulation
library(ggprism)       # For prism-style ggplots
library(ggpubr)        # For geom_signif()


# ----------------------------------------
# Load and preprocess species-level OTU data
# ----------------------------------------

# Set working directory and load OTU data
cat("Loading data ...\n")
setwd("/Users/starklaraliz/Library/CloudStorage/OneDrive-UniversityofCambridge/General - User_Microbiome Function and Diversity/Projects/2024_KlaraStark_HumanMouseGut/scripts/klara/MSvsLC") # <-- Use your actual path
data_otu = read.delim("merged_abundance_table.txt", skip=1L)

# Extract and clean taxonomy information
phylum.df = as.data.frame(data_otu[grepl("s__", data_otu$clade_name), "clade_name"])
colnames(phylum.df) = c("clade_name")
phylum.df = phylum.df[!grepl("t__", phylum.df$clade_name),]
phylum.df = str_split_fixed(phylum.df, ";", 7)
phylum.df = as.data.frame(sub("[dpcofgs]__", "", phylum.df))
colnames(phylum.df) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phylum.df[length(rownames(phylum.df))+1,] = rep("Other", 7)

# Filter OTU table to species-level entries
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = t(species.df[,-1])

# Load and align metadata with species data
metadata = read.delim("MSvsLC_metadata.txt")
rownames(metadata) = metadata$run.accession
metadata.filt = metadata[rownames(species.df), ]

# ----------------------------------------
# Alpha diversity analysis (Shannon-index)
# ----------------------------------------

# Calculate Shannon diversity and merge with metadata
cat("Calculating Shannon-index ...\n")
alpha.div = diversity(species.df, index = "shannon")
metadata.filt$Shannon = alpha.div[rownames(metadata.filt)]
metadata.filt = metadata.filt[metadata.filt$annot1 != "WT", ]
metadata.filt = metadata.filt[metadata.filt$annot1 != "D2 HMA", ]
metadata.filt$annot2 = factor(metadata.filt$annot2, levels = c("luminal content", "mucosal scraping"))

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Shannon ~ sample, data = metadata.filt)
kruskal_p <- kruskal_test_result$p.value

# Post-hoc Dunn’s test and adjust p-values if Kruskal-Wallis is significant
if (kruskal_p < 0.05) {
  dunn_test <- dunn.test(metadata.filt$Shannon, metadata.filt$sample, kw = TRUE, label = TRUE)
  dunn_results <- data.frame(dunn_test$res)
  
  # Adjust p-values using Benjamini-Hochberg correction
  dunn_results$p.adj <- p.adjust(dunn_results$P, method = "BH")
  
  # Create a "Comparison" column by concatenating group names
  dunn_results$Comparison <- paste(dunn_results$Group1, dunn_results$Group2, sep = " - ")
  
  # Filter significant comparisons
  significant_pairs <- dunn_results %>%
    filter(p.adj < 0.05) %>%
    select(Comparison, p.adj)
}

# Prepare list of significant pairs for plotting
signif_pairs <- strsplit(as.character(significant_pairs$comparison), " - ")

# Summarize data for boxplot visualization
summary_df <- metadata.filt %>%
  group_by(annot2) %>%
  summarise(
    median_value = median(Shannon),
    Q1 = quantile(Shannon, 0.25),  # First quartile (Q1)
    Q3 = quantile(Shannon, 0.75)   # Third quartile (Q3)
  )

# ----------------------------------------
# Plotting box plot
# ----------------------------------------

# Generate alpha diversity boxplot with significance annotations
alpha.plot = ggplot(metadata.filt, aes(annot2, Shannon, fill = annot2)) + 
  geom_boxplot(outliers = FALSE, coef = 0, width = 0.5, position = position_dodge(width = 0.5), color = "black") +
  geom_point(position = position_jitter(width = 0.2), colour = "black", size = 2) +
  theme_prism() +
  labs(x = "", y = "Shannon-index") +
  theme(legend.title = element_blank()) +
  theme(axis.title = element_text(size=14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold")) +  # Increase y-axis label size
  theme(axis.text.x = element_text(size=14, face = "bold"),  legend.position = "none") +
  geom_signif(comparisons = list(comp1=c("mucosal scraping", "luminal content")),
              map_signif_level = FALSE, step_increase = 0.1, textsize = 4, y_position = c(3.7),
              color = "black") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # Median as a horizontal line
  scale_fill_manual(values = c("luminal content" = "#FA8", "mucosal scraping" = "#9E2A2F"))  

# Display plot
alpha.plot

