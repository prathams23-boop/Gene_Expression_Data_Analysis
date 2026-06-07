# WGCNA, module preservation analysis
# hub gene identification
# for TCGA-BRCA (Breast Cancer gene expression data)
# dataset cleaned and preprocessed before use

library(WGCNA)
library(tidyverse)
library(gplots)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(dplyr)
library(dendextend)
library(gplots)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(dplyr)
library(GO.db)
library(DESeq2)
library(genefilter)
library(clusterProfiler)

# Load tumor data 

tumor_data <- read.csv("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification/TCGA-BRCA_clean_expression_tumor.csv",
                       row.names=1)
# row.names=1 tells R that the very first column of the 
# spreadsheet contains gene ids and should be used
# to name the rows rather than being treated as 
# numerical data

normal_data <- read.csv("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification/TCGA-BRCA_clean_expression_normal.csv",
                        row.names=1)
dim(tumor_data)
dim(normal_data)

#viewing the first 5 rows and columns
# and seeing column names
tumor_data[1:5, 1:5]
colnames(tumor_data)
normal_data[1:5, 1:5]
colnames(normal_data)

# NORMALISE WITH DESeq2

# 1-create metadata table
metadata_tumor <- data.frame(
  Sample=colnames(tumor_data),
  Condition=rep('Tumor', 1106) # rep is replicate
)

rownames(metadata_tumor) <- metadata_tumor$Sample

metadata_normal <- data.frame(
  Sample=colnames(normal_data),
  Condition=rep('Normal', 113)
)

rownames(metadata_normal) <- metadata_normal$Sample

all(colnames(tumor_data)==metadata_tumor$Sample)
# this MUST return TRUE because column names
# must exactly match metadata row names

all(colnames(normal_data)==metadata_normal$Sample)

# 2-run the normalisation

dds_tumor <- DESeqDataSetFromMatrix(
  countData=round(tumor_data),
  colData=metadata_tumor,
  design=~1
)

dds_tumor <- DESeq(dds_tumor) #running actual DESeq2
normalized_counts_tumor <- counts(dds_tumor, normalized=TRUE)

dds_normal <- DESeqDataSetFromMatrix(
  countData = round(normal_data),
  colData = metadata_normal,
  design = ~1
)

dds_normal <- DESeq(dds_normal)
normalized_counts_normal <- counts(dds_normal, normalized=TRUE)

# applying VST and 
# filtering low-variance genes

vsd_tumor <- varianceStabilizingTransformation(dds_tumor)
vsd_normal <- varianceStabilizingTransformation(dds_normal)

# per gene variance across all samples

rv_tumor <- rowVars(assay(vsd_tumor))
rv_normal <- rowVars(assay(vsd_normal))

# keeping the top 5%

q95_tumor <- quantile(rv_tumor, 0.95)
q95_normal <- quantile(rv_normal, 0.95)

filtered_tumor <- assay(vsd_tumor)[rv_tumor>q95_tumor,]
filtered_normal <- assay(vsd_normal)[rv_normal>q95_normal,]

# aligning (needed when number of tumor is not equal to number of normal)

# 1. Find the genes that survived the 5% filter in BOTH datasets
common_filtered_genes <- intersect(rownames(filtered_tumor), rownames(filtered_normal))

# 2. Subset both filtered matrices to only keep these overlapping genes
filtered_tumor <- filtered_tumor[common_filtered_genes, ]
filtered_normal <- filtered_normal[common_filtered_genes, ]

# 3. Check the final dimensions to make sure they match perfectly
dim(filtered_tumor)
dim(filtered_normal)

# checking the final counts

cat("Genes retained-tumor: ", nrow(filtered_tumor), "\n")
cat("Genes retained-normal: ", nrow(filtered_normal), "\n")

# choosing soft thresholding power

allowWGCNAThreads() # use multiple CPU cores

powers <- c(1:10, seq(from=10, to=20, by=2))

sft_tumor_signed <- pickSoftThreshold(
  t(filtered_tumor), # t is transpose
  powerVector = powers,
  networkType = "signed",
  verbose=5
)

sft_tumor_unsigned <- pickSoftThreshold(
  t(filtered_tumor), # t is transpose
  powerVector = powers,
  networkType = "unsigned",
  verbose=5
)

# the above function previews wat the network
# would look like under every single power 
# level in the vector we defined

sft_tumor_signed
sft_tumor_unsigned
power_tumor_unsigned <- 3
power_tumor_signed_1 <- 7
power_tumor_signed_2 <- 8
power_tumor_signed_3 <- 9

sft_normal_signed <- pickSoftThreshold(
  t(filtered_normal), # t is transpose
  powerVector = powers,
  networkType = "signed",
  verbose=5
)

sft_normal_unsigned <- pickSoftThreshold(
  t(filtered_normal), # t is transpose
  powerVector = powers,
  networkType = "unsigned",
  verbose=5
)

sft_normal_signed
sft_normal_unsigned
power_normal_signed_1 <- 18
power_normal_signed_2 <- 20

power_normal_unsigned_1 <- 9
power_normal_unsigned_2 <- 10


# scale independence and mean connectivity plots


# Set text scaling factor so plot labels are cleanly readable and don't overlap
cex1 <- 0.9



# Split the R graphics device into a 1-row, 2-column window for side-by-side plots
par(mfrow = c(2, 2))

# 1. Tumor Scale Independence Plot (unsigned)
plot(sft_tumor_unsigned$fitIndices[, 1],
     -sign(sft_tumor_unsigned$fitIndices[, 3]) * sft_tumor_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",  # "n" sets up the empty grid boundaries without drawing dots
     main = "Scale independence (Tumor, Unsigned)")

# Populate the empty plot with the actual tested power numbers in red text
text(sft_tumor_unsigned$fitIndices[, 1],
     -sign(sft_tumor_unsigned$fitIndices[, 3]) * sft_tumor_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# Draw the critical horizontal target line at R^2 = 0.90
abline(h = 0.9, col = "red")

# 1. Tumor Scale Independence Plot (signed)
plot(sft_tumor_signed$fitIndices[, 1],
     -sign(sft_tumor_signed$fitIndices[, 3]) * sft_tumor_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",  # "n" sets up the empty grid boundaries without drawing dots
     main = "Scale independence (Tumor,Signed)")

# Populate the empty plot with the actual tested power numbers in red text
text(sft_tumor_signed$fitIndices[, 1],
     -sign(sft_tumor_signed$fitIndices[, 3]) * sft_tumor_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# Draw the critical horizontal target line at R^2 = 0.90
abline(h = 0.9, col = "red")


# 2. Tumor Mean Connectivity Plot (unsigned)
plot(sft_tumor_unsigned$fitIndices[, 1],
     sft_tumor_unsigned$fitIndices[, 5], # Column 5 contains the mean connectivity (mean.k.)
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity (Tumor, Unsigned)")

# Populate the mean connectivity plot with the power numbers
text(sft_tumor_unsigned$fitIndices[, 1],
     sft_tumor_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# 2. Tumor Mean Connectivity Plot (signed)
plot(sft_tumor_signed$fitIndices[, 1],
     sft_tumor_signed$fitIndices[, 5], # Column 5 contains the mean connectivity (mean.k.)
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity (Tumor, Signed)")

# Populate the mean connectivity plot with the power numbers
text(sft_tumor_signed$fitIndices[, 1],
     sft_tumor_signed$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")


# plots for normal tissue

# Reset the graphics layout for the normal tissue plots
par(mfrow = c(2, 2))

# 1. Normal Scale Independence Plot (unsigned)
plot(sft_normal_unsigned$fitIndices[, 1],
     -sign(sft_normal_unsigned$fitIndices[, 3]) * sft_normal_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence (Normal, Unsigned)")
     

# Populate the plot with red numbers
text(sft_normal_unsigned$fitIndices[, 1],
     -sign(sft_normal_unsigned$fitIndices[, 3]) * sft_normal_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# Draw the horizontal target line at R^2 = 0.90
abline(h = 0.9, col = "red")

# 1. Normal Scale Independence Plot (signed)
plot(sft_normal_signed$fitIndices[, 1],
     -sign(sft_normal_signed$fitIndices[, 3]) * sft_normal_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence (Normal, Signed)")

# Populate the plot with red numbers
text(sft_normal_signed$fitIndices[, 1],
     -sign(sft_normal_signed$fitIndices[, 3]) * sft_normal_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# Draw the horizontal target line at R^2 = 0.90
abline(h = 0.9, col = "red")


# 2. Normal Mean Connectivity Plot (unsigned)
plot(sft_normal_unsigned$fitIndices[, 1],
     sft_normal_unsigned$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity (Normal, Unsigned)")

# Populate the mean connectivity plot with numbers
text(sft_normal_unsigned$fitIndices[, 1],
     sft_normal_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# 2. Normal Mean Connectivity Plot (signed)
plot(sft_normal_signed$fitIndices[, 1],
     sft_normal_signed$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity (Normal, Signed)")

# Populate the mean connectivity plot with numbers
text(sft_normal_signed$fitIndices[, 1],
     sft_normal_signed$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")





# CONSTRUCTING CO-EXPRESSION NETWORKS

library(WGCNA)
cor <- WGCNA::cor

#tumor unsigned
brca_netwk_tumor_unsigned <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_unsigned,
  networkType = "unsigned",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "unsigned",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_unsigned",
  verbose=3
  
)

exists("brca_netwk_tumor_unsigned")
saveRDS(brca_netwk_tumor_unsigned, "brca_netwk_tumor_unsigned.rds")

#tumor signed 1
brca_netwk_tumor_signed_1 <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_signed_1,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_signed_1",
  verbose=3
  
)
saveRDS(brca_netwk_tumor_signed_1, "brca_netwk_tumor_signed_1.rds")

#tumor signed 2
brca_netwk_tumor_signed_2 <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_signed_2,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_signed_2",
  verbose=3
  
)
saveRDS(brca_netwk_tumor_signed_2, "brca_netwk_tumor_signed_2.rds")

#tumor signed 3
brca_netwk_tumor_signed_3 <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_signed_3,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_signed_3",
  verbose=3
  
)
saveRDS(brca_netwk_tumor_signed_3, "brca_netwk_tumor_signed_3.rds")

#normal unsigned 1
brca_netwk_normal_unsigned_1 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_unsigned_1,
  networkType = "unsigned",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "unsigned",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_unsigned_1",
  verbose=3
  
)
saveRDS(brca_netwk_normal_unsigned_1, "brca_netwk_normal_unsigned_1.rds")

#normal unsigned 2
brca_netwk_normal_unsigned_2 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_unsigned_2,
  networkType = "unsigned",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "unsigned",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_unsigned_2",
  verbose=3
  
)

saveRDS(brca_netwk_normal_unsigned_2, "brca_netwk_normal_unsigned_2.rds")

#normal signed 1
brca_netwk_normal_signed_1 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_signed_1,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_signed_1",
  verbose=3
  
)

saveRDS(brca_netwk_normal_signed_1, "brca_netwk_normal_signed_1.rds")

#normal signed 2
brca_netwk_normal_signed_2 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_signed_2,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_signed_2",
  verbose=3
  
)

saveRDS(brca_netwk_normal_signed_1, "brca_netwk_normal_signed_1.rds")

# CLUSTER DENDOGRAMS

# convert raw wgcna module names to colors
colors_tum_signed_1   <- labels2colors(brca_netwk_tumor_signed_1$colors)
colors_tum_signed_2 <- labels2colors(brca_netwk_tumor_signed_2$colors)
colors_tum_signed_3 <- labels2colors(brca_netwk_tumor_signed_3$colors)
colors_tum_unsigned <- labels2colors(brca_netwk_tumor_unsigned$colors)

colors_norm_signed_1  <- labels2colors(brca_netwk_normal_signed_1$colors)
colors_norm_signed_2  <- labels2colors(brca_netwk_normal_signed_2$colors)
colors_norm_unsigned_1 <- labels2colors(brca_netwk_normal_unsigned_1$colors)
colors_norm_unsigned_2 <- labels2colors(brca_netwk_normal_unsigned_2$colors)


# plot tumor networks signed vs unsigned side by side
# Set up a 2-row, 2-column plotting grid
par(mfrow = c(2, 2))

#  Tumor Signed Dendrogram 1
plotDendroAndColors(
  brca_netwk_tumor_signed_1$dendrograms[[1]],                  # Raw clustering tree structures
  colors_tum_signed_1[brca_netwk_tumor_signed_1$blockGenes[[1]]], # Aligns colors to match correct leaves
  "Module colors",                                      # Title of the color track bar
  dendroLabels = FALSE,                                 # Prevents overlapping gene names text
  hang = 0.03,                                          # Clean visual alignment above color bar
  addGuide = TRUE,                                      # Adds vertical guide lines
  guideHang = 0.05,
  main = "Tumor Network - SIGNED_1"
)

#  Tumor Signed Dendrogram 2
plotDendroAndColors(
  brca_netwk_tumor_signed_2$dendrograms[[1]],                  # Raw clustering tree structures
  colors_tum_signed_2[brca_netwk_tumor_signed_2$blockGenes[[1]]], # Aligns colors to match correct leaves
  "Module colors",                                      # Title of the color track bar
  dendroLabels = FALSE,                                 # Prevents overlapping gene names text
  hang = 0.03,                                          # Clean visual alignment above color bar
  addGuide = TRUE,                                      # Adds vertical guide lines
  guideHang = 0.05,
  main = "Tumor Network - SIGNED_2"
)

#  Tumor Signed Dendrogram 3
plotDendroAndColors(
  brca_netwk_tumor_signed_3$dendrograms[[1]],                  # Raw clustering tree structures
  colors_tum_signed_3[brca_netwk_tumor_signed_3$blockGenes[[1]]], # Aligns colors to match correct leaves
  "Module colors",                                      # Title of the color track bar
  dendroLabels = FALSE,                                 # Prevents overlapping gene names text
  hang = 0.03,                                          # Clean visual alignment above color bar
  addGuide = TRUE,                                      # Adds vertical guide lines
  guideHang = 0.05,
  main = "Tumor Network - SIGNED_3"
)

# Tumor Unsigned Dendrogram 
plotDendroAndColors(
  brca_netwk_tumor_unsigned$dendrograms[[1]],
  colors_tum_unsigned[brca_netwk_tumor_unsigned$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Tumor Network - UNSIGNED"
)


# same plots for normal


# normal signed 1
plotDendroAndColors(
  brca_netwk_normal_signed_1$dendrograms[[1]],
  colors_norm_signed_1[brca_netwk_normal_signed_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - SIGNED_1"
)

# normal signed 2
plotDendroAndColors(
  brca_netwk_normal_signed_2$dendrograms[[1]],
  colors_norm_signed_2[brca_netwk_normal_signed_2$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - SIGNED_2"
)




# normal unsigned 1
plotDendroAndColors(
  brca_netwk_normal_unsigned_1$dendrograms[[1]],
  colors_norm_unsigned_1[brca_netwk_normal_unsigned_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - UNSIGNED_1"
)

# normal unsigned 2
plotDendroAndColors(
  brca_netwk_normal_unsigned_2$dendrograms[[1]],
  colors_norm_unsigned_2[brca_netwk_normal_unsigned_2$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - UNSIGNED_2"
)

# Reset plotting grid back to default
par(mfrow = c(1, 1))

# **IMPORTANT** 
# now based on these dendrograms, we figure out which of the candidate powers
# is actually the most suitable and then we go ahead with only 1 network for 
# each signed and unsigned tumor and normal for the further analysis


# MODULE SIZE TABLES

table(netwk_tumor_unsigned$colors)
table(netwk_tumor_signed$colors)
table(netwk_normal_unsigned$colors)
table(netwk_normal_signed$colors)

# MODULE EIGENGENE COMPARISION
# signed and unsigned


# ------------------------------------------------------------------------------
# STEP 1: Calculate and Clean Eigengenes for All 4 Networks
# ------------------------------------------------------------------------------

# --- 1. Tumor Signed
MEs_tum_signed <- moduleEigengenes(t(filtered_tumor), colors = colors_tum_signed)$eigengenes
MEs_tum_signed <- orderMEs(MEs_tum_signed)
colnames(MEs_tum_signed) <- gsub("ME", "", colnames(MEs_tum_signed))

# --- 2. Tumor Unsigned
MEs_tum_unsigned <- moduleEigengenes(t(filtered_tumor), colors = colors_tum_unsigned)$eigengenes
MEs_tum_unsigned <- orderMEs(MEs_tum_unsigned)
colnames(MEs_tum_unsigned) <- gsub("ME", "", colnames(MEs_tum_unsigned))

# --- 3. Normal Signed
MEs_norm_signed <- moduleEigengenes(t(filtered_normal), colors = colors_norm_signed)$eigengenes
MEs_norm_signed <- orderMEs(MEs_norm_signed)
colnames(MEs_norm_signed) <- gsub("ME", "", colnames(MEs_norm_signed))

# --- 4. Normal Unsigned
MEs_norm_unsigned <- moduleEigengenes(t(filtered_normal), colors = colors_norm_unsigned)$eigengenes
MEs_norm_unsigned <- orderMEs(MEs_norm_unsigned)
colnames(MEs_norm_unsigned) <- gsub("ME", "", colnames(MEs_norm_unsigned))


# ------------------------------------------------------------------------------
# STEP 2: Compute Pairwise Pearson Correlation Matrices
# ------------------------------------------------------------------------------
cor_tum_signed    <- cor(MEs_tum_signed)
cor_tum_unsigned  <- cor(MEs_tum_unsigned)
cor_norm_signed   <- cor(MEs_norm_signed)
cor_norm_unsigned <- cor(MEs_norm_unsigned)


# ------------------------------------------------------------------------------
# STEP 3: Verify Eigengene Dimensions
# ------------------------------------------------------------------------------
cat("=== Tumor Signed Matrix Dimensions ===\n")
print(dim(cor_tum_signed))

cat("\n=== Tumor Unsigned Matrix Dimensions ===\n")
print(dim(cor_tum_unsigned))


# ==============================================================================
# LISTING 21: GENERATE FIGURE 9 — EIGENGENE CORRELATION HEATMAPS
# ==============================================================================

# Ensure the gplots library is loaded for heatmap.2
library(gplots)

# Step 1: Calculate the 4 internal Pearson correlation matrices
cor_tum_unsigned  <- cor(MEs_tum_unsigned)
cor_tum_signed    <- cor(MEs_tum_signed)
cor_norm_unsigned <- cor(MEs_norm_unsigned)
cor_norm_signed   <- cor(MEs_norm_signed)

# Step 2: Split the plotting window into a 2x2 grid to see all 4 at once
par(mfrow = c(2, 2))

# --- HEATMAP 9A: Tumor Unsigned ---
heatmap.2(cor_tum_unsigned,
          main = "Tumor Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# --- HEATMAP 9B: Tumor Signed ---
heatmap.2(cor_tum_signed,
          main = "Tumor Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# --- HEATMAP 9C: Normal Unsigned ---
heatmap.2(cor_norm_unsigned,
          main = "Normal Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# --- HEATMAP 9D: Normal Signed ---
heatmap.2(cor_norm_signed,
          main = "Normal Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# Reset plotting grid back to default single panel
par(mfrow = c(1, 1))


# ==============================================================================
# SECTION 11: EXPORT ALL 4 NETWORKS TO CYTOSCAPE FORMAT
# ==============================================================================
# 1. Install BiocManager from CRAN (if you don't have it already)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Use BiocManager to download and install the genome database
BiocManager::install("org.Hs.eg.db")
# Ensure required database libraries are active for gene symbol translation
library(org.Hs.eg.db)
library(tidyverse)

# Create a clean folder in your directory to store the output text files
dir.create("modules", showWarnings = FALSE) [cite: 18]

# ------------------------------------------------------------------------------
# FUNCTION: Reusable pipeline to extract Edge and Node files
# ------------------------------------------------------------------------------
export_to_cytoscape <- rm_duplicate_and_save <- function(filtered_data, power_val, net_type, colors_vector, prefix) {
  
  cat("...Processing TOM for", prefix, "...\n")
  # 1. Calculate the Topological Overlap Matrix
  TOM_mat <- TOMsimilarityFromExpr(t(filtered_data), power = power_val, networkType = net_type)
  row.names(TOM_mat) <- row.names(filtered_data)
  colnames(TOM_mat) <- row.names(filtered_data)
  
  # 2. Reshape into long-format Edge List and filter weak connections (> 0.1)
  edge_list <- data.frame(TOM_mat) %>%
    mutate(gene1 = row.names(.)) %>%
    pivot_longer(-gene1, names_to = "gene2", values_to = "weight") %>%
    filter(gene1 != gene2, weight > 0.1)
  
  # 3. Remove duplicate undirected tracking pairs (A-B vs B-A)
  edge_list <- edge_list %>%
    mutate(pair = pmap_chr(list(gene1, gene2), ~ paste(sort(c(...)), collapse = "_"))) %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  # 4. Translate Ensembl IDs to official human Gene Symbols for both columns
  cat("...Translating Gene Symbols for", prefix, "...\n")
  edge_list$gene1.name <- mapIds(org.Hs.eg.db, keys = edge_list$gene1, column = "SYMBOL", keytype = "ENSEMBL")
  edge_list$gene2.name <- mapIds(org.Hs.eg.db, keys = edge_list$gene2, column = "SYMBOL", keytype = "ENSEMBL")
  
  # Fill in missing symbols with original Ensembl IDs as a fallback safety check
  edge_list$gene1.name <- ifelse(is.na(edge_list$gene1.name), edge_list$gene1, edge_list$gene1.name)
  edge_list$gene2.name <- ifelse(is.na(edge_list$gene2.name), edge_list$gene2, edge_list$gene2.name)
  
  # 5. Construct the corresponding Node List mapping genes to their final colors
  node_list <- data.frame(
    GeneID = row.names(filtered_data),
    ModuleColor = colors_vector
  )
  node_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = node_list$GeneID, column = "SYMBOL", keytype = "ENSEMBL")
  node_list$GeneSymbol <- ifelse(is.na(node_list$GeneSymbol), node_list$GeneID, node_list$GeneSymbol)
  
  # 6. Save data tables to your "modules" folder as clean, tab-separated text files
  write.table(edge_list, file = paste0("modules/Edges_", prefix, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(node_list, file = paste0("modules/Nodes_", prefix, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("=== Successfully exported Edge and Node files for:", prefix, "===\n\n")
}

# ------------------------------------------------------------------------------
# EXECUTE RUNS FOR ALL 4 INDIVIDUAL CONFIGURATIONS
# ------------------------------------------------------------------------------

# Run 1: Tumor Unsigned
export_to_cytoscape(filtered_tumor, power_tumor_unsigned, "unsigned", colors_tum_unsigned, "Tumor_Unsigned")

# Run 2: Tumor Signed
export_to_cytoscape(filtered_tumor, power_tumor_signed, "signed", colors_tum_signed, "Tumor_Signed")

# Run 3: Normal Unsigned
export_to_cytoscape(filtered_normal, power_normal_unsigned, "unsigned", colors_norm_unsigned, "Normal_Unsigned")

# Run 4: Normal Signed
export_to_cytoscape(filtered_normal, power_normal_signed, "signed", colors_norm_signed, "Normal_Signed")


# ==============================================================================
# SECTION 12: MODULE PRESERVATION ANALYSIS (TUMOR AS REFERENCE)
# ==============================================================================
library(WGCNA)
library(tidyverse)
allowWGCNAThreads() # Accelerates calculation using multiple CPU cores

# ------------------------------------------------------------------------------
# BLOCK 1: SIGNED NETWORK MODULE PRESERVATION
# ------------------------------------------------------------------------------
cat("\n=== STEP 1: Preparing Inputs for SIGNED Network ===\n")

# 1. Package the Transposed Data Matrices (Samples as Rows, Genes as Columns)
multiData_signed <- list(
  Tumor  = list(data = t(filtered_tumor)),
  Normal = list(data = t(filtered_normal))
)

# 2. Extract the SIGNED Tumor & Normal Character Color Vectors
tumor_colors_signed  <- labels2colors(netwk_tumor_signed$colors)
names(tumor_colors_signed) <- names(netwk_tumor_signed$colors)

normal_colors_signed <- labels2colors(netwk_normal_signed$colors)
names(normal_colors_signed) <- names(netwk_normal_signed$colors)

multiColor_signed <- list(
  Tumor  = tumor_colors_signed,
  Normal = normal_colors_signed
)

# 3. Quick Integrity Check (Both must return TRUE)
if(all(names(multiColor_signed$Tumor) %in% rownames(filtered_tumor)) && 
   all(names(multiColor_signed$Normal) %in% rownames(filtered_normal))) {
  cat("✓ Sanity checks passed for Signed Network. Proceeding to calculations...\n")
} else {
  stop("✕ Error: Gene name mismatch detected in Signed Network structures!")
}

# 4. Run the Permutation Math for Signed
cat("\n...Running SIGNED Module Preservation (100 loops)...\n")
cat("Note: This is heavily CPU-bound and will take a significant amount of time.\n")

preservation_signed <- modulePreservation(
  multiData         = multiData_signed,
  multiColor        = multiColor_signed,
  referenceNetworks = 1,       # 1 = Tumor is the Reference framework
  nPermutations     = 100,     # 100 baseline loops for statistical rigor
  randomSeed        = 12345,   # Fixes the seed for identical reproducibility
  verbose           = 3
)

# 5. Extract Statistics and Filter out the Artificial "Gold" Control Module
stats_signed <- preservation_signed$preservation$Z$ref.Tumor$inColumnsAlsoPresentIn.Normal
stats_signed <- stats_signed[rownames(stats_signed) != "gold", ]

cat("\n=== SIGNED Network Permutations Complete! ===\n")
print(head(stats_signed))


# ------------------------------------------------------------------------------
# BLOCK 2: UNSIGNED NETWORK MODULE PRESERVATION
# ------------------------------------------------------------------------------
cat("\n=== STEP 2: Preparing Inputs for UNSIGNED Network ===\n")

# 1. Package Data Matrices
multiData_unsigned <- list(
  Tumor  = list(data = t(filtered_tumor)),
  Normal = list(data = t(filtered_normal))
)

# 2. Extract the UNSIGNED Tumor & Normal Character Color Vectors
tumor_colors_unsigned  <- labels2colors(netwk_tumor_unsigned$colors)
names(tumor_colors_unsigned) <- names(netwk_tumor_unsigned$colors)

normal_colors_unsigned <- labels2colors(netwk_normal_unsigned$colors)
names(normal_colors_unsigned) <- names(netwk_normal_unsigned$colors)

multiColor_unsigned <- list(
  Tumor  = tumor_colors_unsigned,
  Normal = normal_colors_unsigned
)

# 3. Quick Integrity Check (Both must return TRUE)
if(all(names(multiColor_unsigned$Tumor) %in% rownames(filtered_tumor)) && 
   all(names(multiColor_unsigned$Normal) %in% rownames(filtered_normal))) {
  cat("✓ Sanity checks passed for Unsigned Network. Proceeding to calculations...\n")
} else {
  stop("✕ Error: Gene name mismatch detected in Unsigned Network structures!")
}

# 4. Run the Permutation Math for Unsigned
cat("\n...Running UNSIGNED Module Preservation (100 loops)...\n")

preservation_unsigned <- modulePreservation(
  multiData         = multiData_unsigned,
  multiColor        = multiColor_unsigned,
  referenceNetworks = 1,       # 1 = Tumor is the Reference framework
  nPermutations     = 100,     
  randomSeed        = 12345,   
  verbose           = 3
)

# 5. Extract Statistics and Filter out the Artificial "Gold" Control Module
stats_unsigned <- preservation_unsigned$preservation$Z$ref.Tumor$inColumnsAlsoPresentIn.Normal
stats_unsigned <- stats_unsigned[rownames(stats_unsigned) != "gold", ]

cat("\n=== UNSIGNED Network Permutations Complete! ===\n")
print(head(stats_unsigned))

# ==============================================================================
# SECTION 12.4 & 12.5: EXTRACT STATS & GENERATE Z-SUMMARY DOT PLOTS
# ==============================================================================
library(ggplot2)

# ------------------------------------------------------------------------------
# PART A: Plotting the SIGNED Network Preservation
# ------------------------------------------------------------------------------
# Create a clean dataframe for ggplot
plot_data_signed <- data.frame(
  Module    = rownames(stats_signed),
  Z_summary = as.numeric(stats_signed$Zsummary.pres)
)

# Generate the Figure 20 Dot Plot for Signed
ggplot(plot_data_signed, aes(x = Z_summary, y = reorder(Module, Z_summary))) +
  geom_point(aes(color = Module), size = 5) +
  scale_color_identity() + # Uses the literal module text color for the dots
  geom_vline(xintercept = 2, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Module Preservation Z-summary (SIGNED Network)",
    subtitle = "Reference: Tumor | Test: Normal",
    x = "Z-summary Score",
    y = "Tumor Modules"
  ) +
  theme_minimal() +
  annotate("text", x = 1, y = 1, label = "No Preservation", color = "blue", angle = 90, vjust = -0.5) +
  annotate("text", x = 6, y = 1, label = "Moderate", color = "purple", angle = 90, vjust = -0.5) +
  annotate("text", x = 12, y = 1, label = "Strong Preservation", color = "red", angle = 90, vjust = -0.5)

# ------------------------------------------------------------------------------
# PART B: Plotting the UNSIGNED Network Preservation
# ------------------------------------------------------------------------------
# Create a clean dataframe for ggplot
plot_data_unsigned <- data.frame(
  Module    = rownames(stats_unsigned),
  Z_summary = as.numeric(stats_unsigned$Zsummary.pres)
)

# Generate the Figure 20 Dot Plot for Unsigned
ggplot(plot_data_unsigned, aes(x = Z_summary, y = reorder(Module, Z_summary))) +
  geom_point(aes(color = Module), size = 5) +
  scale_color_identity() + 
  geom_vline(xintercept = 2, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Module Preservation Z-summary (UNSIGNED Network)",
    subtitle = "Reference: Tumor | Test: Normal",
    x = "Z-summary Score",
    y = "Tumor Modules"
  ) +
  theme_minimal() +
  annotate("text", x = 1, y = 1, label = "No Preservation", color = "blue", angle = 90, vjust = -0.5) +
  annotate("text", x = 6, y = 1, label = "Moderate", color = "purple", angle = 90, vjust = -0.5) +
  annotate("text", x = 12, y = 1, label = "Strong Preservation", color = "red", angle = 90, vjust = -0.5)


# ==============================================================================
# SECTION 12.6: EXTRACT GENE LISTS PER MODULE (FIGURE 21)
# ==============================================================================

# 1. Get all unique module color names from your Signed network construction
module_names_signed <- unique(tumor_colors_signed)

# 2. Group the Ensembl Gene IDs into named list buckets based on their color
genes_in_modules_signed <- lapply(module_names_signed, function(module) {
  names(tumor_colors_signed[tumor_colors_signed == module])
})

# 3. Attach the official color names to each list bucket
names(genes_in_modules_signed) <- module_names_signed

# 4. Print the structure to your console (Generates Figure 21 Output)
cat("\n=== EXTRACTED MODULE GENE COUNTS (SIGNED WORKFLOW) ===\n")
str(genes_in_modules_signed)

# ==============================================================================
# SECTION 13: GO ENRICHMENT ANALYSIS (FIGURES 22-25)
# ==============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------------------------
# STEP 1: Set Up Clean Output Directories (Listing 31)
# ------------------------------------------------------------------------------
dir.create("enrich", showWarnings = FALSE)
dir.create("enrich/GO_T_N", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# STEP 2: Build a Reusable GO Enrichment Function (Listing 32)
# ------------------------------------------------------------------------------
perform_go_enrichment <- function(gene_list, ontology, output_path) {
  
  cat(paste("\n...Running GO Enrichment for Ontology:", ontology, "...\n"))
  
  go_results <- enrichGO(
    gene          = gene_list,      # Our vector of Ensembl Gene IDs from the module
    OrgDb         = org.Hs.eg.db,   # The human genome annotation database mapping file
    keyType       = "ENSEMBL",      # Explicitly states our input format is Ensembl IDs
    ont           = ontology,       # Can accept "BP", "CC", or "MF"
    pAdjustMethod = "BH",           # Benjamini-Hochberg false-discovery rate correction
    pvalueCutoff  = 0.05,           # Statistical significance threshold
    qvalueCutoff  = 0.05            # Strict false-positive filtering cutoff
  )
  
  # Save raw results automatically to a spreadsheet for your report appendices
  write.csv(as.data.frame(go_results), file = output_path, row.names = TRUE)
  
  return(go_results)
}

# ------------------------------------------------------------------------------
# STEP 3: Run GO Analysis for the Blue Module (Listing 33)
# ------------------------------------------------------------------------------
# 1. Isolate your target blue module genes (ensuring background 'grey' isn't mixed in)
blue_genes <- genes_in_modules_signed$blue

# 2. Compute Biological Processes (BP) for Figure 22
go_blue_bp <- perform_go_enrichment(
  gene_list   = blue_genes, 
  ontology    = "BP", 
  output_path = "enrich/GO_T_N/GO_BP_blue.csv"
)

# 3. Compute Cellular Components (CC)
go_blue_cc <- perform_go_enrichment(
  gene_list   = blue_genes, 
  ontology    = "CC", 
  output_path = "enrich/GO_T_N/GO_CC_blue.csv"
)

# 4. Compute Molecular Functions (MF)
go_blue_mf <- perform_go_enrichment(
  gene_list   = blue_genes, 
  ontology    = "MF", 
  output_path = "enrich/GO_T_N/GO_MF_blue.csv"
)

# ------------------------------------------------------------------------------
# STEP 4: Visualize the Top Biological Processes (Figure 22 Output)
# ------------------------------------------------------------------------------
cat("\n=== Generating Figure 22: Blue Module Enrichment Plot ===\n")
dotplot(go_blue_bp, showCategory = 10, title = "Top 10 Enriched Biological Processes (Blue Module)")


# ==============================================================================
# SECTION 13.2 (CONTINUED): AUTOMATED GO ENRICHMENT FOR ALL MODULES (LISTING 34)
# ==============================================================================

# 1. Gather all unique module names and strip out the "grey" noise module
module_names_to_test <- names(genes_in_modules_signed)
module_names_to_test <- module_names_to_test[module_names_to_test != "grey"]

cat("Found modules to analyze:", paste(module_names_to_test, collapse=", "), "\n")

# 2. Run the loop engine across all remaining colors
for (mod in module_names_to_test) {
  
  # Extract the specific gene list vector for the current loop color
  gene_list <- genes_in_modules_signed[[mod]]
  
  # Cycle through all 3 branches of the Gene Ontology dictionary
  for (ont in c("BP", "CC", "MF")) {
    
    # Dynamically build the file path (e.g., "enrich/GO_T_N/GO_BP_turquoise.csv")
    out_path <- paste0("enrich/GO_T_N/GO_", ont, "_", mod, ".csv")
    
    message("Running GO-", ont, " for the ", mod, " module...")
    
    # Use tryCatch so if a tiny module has zero significant terms, the loop doesn't crash
    tryCatch({
      perform_go_enrichment(
        gene_list   = gene_list, 
        ontology    = ont, 
        output_path = out_path
      )
    }, error = function(e) {
      message("  ℹ Note: No statistically significant terms found for: ", mod, " (", ont, ")")
    })
  }
}

cat("\n=== All Module GO Enrichments Computed and Saved to 'enrich/GO_T_N/'! ===\n")


# ==============================================================================
# SECTION 13.3 (ALTERNATIVE): INDIVIDUAL DOTPLOTS FOR EVERY MODULE
# ==============================================================================
library(clusterProfiler)
library(ggplot2)

# 1. Gather all unique module names and strip out the "grey" noise
module_names_to_plot <- names(genes_in_modules_signed)
module_names_to_plot <- module_names_to_plot[module_names_to_plot != "grey"]

# 2. Loop through each color and generate its unique individual plot
for (mod in module_names_to_plot) {
  
  # Dynamically locate the saved CSV file we generated in the previous step
  csv_file <- paste0("enrich/GO_T_N/GO_BP_", mod, ".csv")
  
  # Check if the file exists and has data before plotting
  if (file.exists(csv_file)) {
    go_data <- read.csv(csv_file)
    
    if (nrow(go_data) > 0) {
      cat("Generating individual Biological Process plot for module:", mod, "\n")
      
      # Read the data back into a formal enrichResult object so clusterProfiler can plot it
      # If your original go objects are still in your R environment, we can plot them directly:
      obj_name <- paste0("go_", mod, "_bp")
      
      # Open a new plotting window so they don't overwrite each other in RStudio
      if (.Platform$OS.type == "windows") { dev.new() } else { x11() }
      
      # Construct the plot title dynamically based on the current color
      plot_title <- paste0("Top 10 Enriched Biological Processes (", toupper(mod), " Module)")
      
      # Use tryCatch in case a module has too few terms to render safely
      tryCatch({
        # We read the file data directly into a ggplot-friendly format for custom standalone plots
        # Sorting by significance (p.adjust) and grabbing the top 10 rows
        top_10 <- head(go_data[order(go_data$p.adjust), ], 10)
        
        # Calculate GeneRatio as a decimal number for the X-axis mapping
        # Splitting strings like "5/300" into a numeric fraction
        top_10$GeneRatio_num <- sapply(top_10$GeneRatio, function(x) {
          num_den <- as.numeric(strsplit(x, "/")[[1]])
          return(num_den[1] / num_den[2])
        })
        
        # Generate the crisp standalone dotplot
        p <- ggplot(top_10, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue") +
          labs(
            title = plot_title,
            subtitle = "Reference: Tumor Network Construction",
            x = "Gene Ratio",
            y = "Biological Process",
            color = "Adjusted p-value",
            size = "Gene Count"
          ) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        
        print(p)
        
        # Automatically save each image as a crisp PNG file inside your enrich folder!
        ggsave(filename = paste0("enrich/GO_T_N/dotplot_standalone_", mod, ".png"), 
               plot = p, width = 8, height = 6, dpi = 300)
        
      }, error = function(e) {
        message("  ℹ Could not plot standalone figure for: ", mod)
      })
    }
  }
}
