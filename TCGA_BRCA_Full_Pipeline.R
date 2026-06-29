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

tumor_data <- read.csv("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run02/TCGA-BRCA_protein_coding_tumor.csv",
                       row.names = 1,
                       check.names = FALSE)

normal_data <- read.csv("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run02/TCGA-BRCA_protein_coding_normal.csv",
                        row.names = 1,
                        check.names = FALSE)

# The first row of these CSVs is gene_name, which was column 2 in the
# extraction script. After row.names=1 takes ensembl_id as rownames,
# gene_name becomes the first data column. We need to drop it before
# passing numeric data to DESeq2.
# Save the gene name mapping first, then drop that column.

gene_name_map <- data.frame(
  ensembl_id = rownames(tumor_data),
  gene_name  = tumor_data[, 1],
  stringsAsFactors = FALSE
)

tumor_data  <- tumor_data[, -1]   # drop gene_name column
normal_data <- normal_data[, -1]  # drop gene_name column

dim(tumor_data)
dim(normal_data)

# viewing the first 5 rows and columns
# and seeing column names
tumor_data[1:5, 1:5]
colnames(tumor_data)[1:5]
normal_data[1:5, 1:5]
colnames(normal_data)[1:5]


# differential expression analysis

# 1 - build combined matrix and metadata
combined_counts <- cbind(tumor_data, normal_data)

metadata_combined <- data.frame(
  Sample    = c(colnames(tumor_data), colnames(normal_data)),
  Condition = c(rep("Tumor",  ncol(tumor_data)),
                rep("Normal", ncol(normal_data)))
)
rownames(metadata_combined) <- metadata_combined$Sample
metadata_combined$Condition  <- factor(metadata_combined$Condition,
                                       levels = c("Normal", "Tumor"))

# safety check
all(colnames(combined_counts) == metadata_combined$Sample)
# must return TRUE

# 2 - build DESeqDataSet with full design
dds_combined <- DESeqDataSetFromMatrix(
  countData = round(combined_counts),
  colData   = metadata_combined,
  design    = ~Condition
)

# 3 - run DE analysis
dds_combined <- DESeq(dds_combined)

# 4 - extract results: Tumor vs Normal
# lfcThreshold=1 means we require at least 2-fold change
de_results <- results(dds_combined,
                      contrast      = c("Condition", "Tumor", "Normal"),
                      alpha         = 0.05)

summary(de_results)

# 5 - filter to significant DEGs
# padj < 0.05 and absolute log2FC >= 1
de_results_df <- as.data.frame(de_results)
de_results_df <- de_results_df[!is.na(de_results_df$padj), ]

sig_degs <- rownames(de_results_df[
  de_results_df$padj < 0.05 & abs(de_results_df$log2FoldChange) >= 1, 
])

cat("Number of significant DEGs:", length(sig_degs), "\n")

# save the full DE results table for reference
write.csv(de_results_df, "TCGA-BRCA_DE_results_tumor_vs_normal.csv")
write.csv(data.frame(ensembl_id = sig_degs), "TCGA-BRCA_sig_DEGs_list.csv",
          row.names = FALSE)


# NORMALISE WITH DESeq2 (separately, design=~1, for downstream VST)
# We normalise each tissue type independently for WGCNA, using only the DEG set


# 1 - create metadata tables for each tissue separately
metadata_tumor <- data.frame(
  Sample    = colnames(tumor_data),
  Condition = rep("Tumor", ncol(tumor_data))
)
rownames(metadata_tumor) <- metadata_tumor$Sample

metadata_normal <- data.frame(
  Sample    = colnames(normal_data),
  Condition = rep("Normal", ncol(normal_data))
)
rownames(metadata_normal) <- metadata_normal$Sample

all(colnames(tumor_data)  == metadata_tumor$Sample)
# must return TRUE
all(colnames(normal_data) == metadata_normal$Sample)
# must return TRUE

# 2 - subset both matrices to DEGs only before building the DESeq objects
tumor_deg  <- tumor_data[rownames(tumor_data)   %in% sig_degs, ]
normal_deg <- normal_data[rownames(normal_data) %in% sig_degs, ]

cat("DEGs present in tumor matrix: ",  nrow(tumor_deg),  "\n")
cat("DEGs present in normal matrix: ", nrow(normal_deg), "\n")

# 3 - run normalisation (design=~1: no group comparison, just size factor estimation)
dds_tumor <- DESeqDataSetFromMatrix(
  countData = round(tumor_deg),
  colData   = metadata_tumor,
  design    = ~1
)
dds_tumor <- DESeq(dds_tumor)
normalized_counts_tumor <- counts(dds_tumor, normalized = TRUE)

dds_normal <- DESeqDataSetFromMatrix(
  countData = round(normal_deg),
  colData   = metadata_normal,
  design    = ~1
)
dds_normal <- DESeq(dds_normal)
normalized_counts_normal <- counts(dds_normal, normalized = TRUE)


# applying VST and 
# filtering low-variance genes

vsd_tumor <- varianceStabilizingTransformation(dds_tumor)
vsd_normal <- varianceStabilizingTransformation(dds_normal)

# the DEG set is already biologically filtered
# just apply VST and use all DEGs

filtered_tumor  <- assay(vsd_tumor)
filtered_normal <- assay(vsd_normal)

# align to common genes (should already be identical but check anyway)
common_filtered_genes <- intersect(rownames(filtered_tumor), rownames(filtered_normal))
filtered_tumor  <- filtered_tumor[common_filtered_genes, ]
filtered_normal <- filtered_normal[common_filtered_genes, ]

dim(filtered_tumor)
dim(filtered_normal)
cat("Genes going into WGCNA:", nrow(filtered_tumor), "\n")

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


power_tumor_unsigned <- 4
power_tumor_signed <- 9


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

power_normal_signed_1 <- 14 # finalised
power_normal_signed_2 <- 16

power_normal_unsigned_1 <- 7 # finalised
power_normal_unsigned_2 <- 8


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


# since module size was bit high and number of genes were lesser in each module
# we will do some hyperparameter tuning

brca_netwk_tumor_unsigned_1 <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_unsigned,
  networkType = "unsigned",
  deepSplit = 1,
  minModuleSize = 50,
  TOMType = "unsigned",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_unsigned_1",
  verbose=3
  
)



#tumor signed 
brca_netwk_tumor_signed <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_signed,
  networkType = "signed",
  deepSplit = 2,
  minModuleSize = 30,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_signed",
  verbose=3
  
)

# same tuning we will try here

brca_netwk_tumor_signed_1 <- blockwiseModules(
  t(filtered_tumor),
  power=power_tumor_signed,
  networkType = "signed",
  deepSplit = 1,
  minModuleSize = 50,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "brca_tumor_signed_1",
  verbose=3
  
)



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

# tuning

#normal unsigned 1.1
brca_netwk_normal_unsigned_1_1 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_unsigned_1,
  networkType = "unsigned",
  deepSplit = 1,
  minModuleSize = 50,
  TOMType = "unsigned",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_unsigned_1_1",
  verbose=3
  
)

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

# tuning 

#normal signed 1.1
brca_netwk_normal_signed_1_1 <- blockwiseModules(
  t(filtered_normal),
  power=power_normal_signed_1,
  networkType = "signed",
  deepSplit = 1,
  minModuleSize = 50,
  TOMType = "signed",
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "normal_signed_1_1",
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

# MODULE SIZE TABLES

table(brca_netwk_tumor_unsigned$colors)
table(brca_netwk_tumor_unsigned_1$colors)
table(brca_netwk_tumor_signed$colors)
table(brca_netwk_tumor_signed_1$colors)



table(brca_netwk_normal_signed_1$colors)
table(brca_netwk_normal_signed_1_1$colors)

table(brca_netwk_normal_signed_2$colors)

table(brca_netwk_normal_unsigned_1$colors)
table(brca_netwk_normal_unsigned_1_1$colors)

table(brca_netwk_normal_unsigned_2$colors)

# CLUSTER DENDOGRAMS

# convert raw wgcna module names to colors
colors_tum_signed   <- labels2colors(brca_netwk_tumor_signed$colors)
colors_tum_signed_1 <- labels2colors(brca_netwk_tumor_signed_1$colors)

colors_tum_unsigned <- labels2colors(brca_netwk_tumor_unsigned$colors)
colors_tum_unsigned_1 <- labels2colors(brca_netwk_tumor_unsigned_1$colors)

colors_norm_signed_1  <- labels2colors(brca_netwk_normal_signed_1$colors)
colors_norm_signed_1_1  <- labels2colors(brca_netwk_normal_signed_1_1$colors)
colors_norm_signed_2  <- labels2colors(brca_netwk_normal_signed_2$colors)
colors_norm_unsigned_1 <- labels2colors(brca_netwk_normal_unsigned_1$colors)
colors_norm_unsigned_1_1 <- labels2colors(brca_netwk_normal_unsigned_1_1$colors)
colors_norm_unsigned_2 <- labels2colors(brca_netwk_normal_unsigned_2$colors)


# plot tumor networks signed vs unsigned side by side
# Set up a 2-row, 2-column plotting grid
par(mfrow = c(2, 2))

#  Tumor Signed Dendrogram 
plotDendroAndColors(
  brca_netwk_tumor_signed$dendrograms[[1]],                  # Raw clustering tree structures
  colors_tum_signed[brca_netwk_tumor_signed$blockGenes[[1]]], # Aligns colors to match correct leaves
  "Module colors",                                      # Title of the color track bar
  dendroLabels = FALSE,                                 # Prevents overlapping gene names text
  hang = 0.03,                                          # Clean visual alignment above color bar
  addGuide = TRUE,                                      # Adds vertical guide lines
  guideHang = 0.05,
  main = "Tumor Network - SIGNED_1"
)

#  Tumor Signed Dendrogram (tuned)
plotDendroAndColors(
  brca_netwk_tumor_signed_1$dendrograms[[1]],                  # Raw clustering tree structures
  colors_tum_signed_1[brca_netwk_tumor_signed_1$blockGenes[[1]]], # Aligns colors to match correct leaves
  "Module colors",                                      # Title of the color track bar
  dendroLabels = FALSE,                                 # Prevents overlapping gene names text
  hang = 0.03,                                          # Clean visual alignment above color bar
  addGuide = TRUE,                                      # Adds vertical guide lines
  guideHang = 0.05,
  main = "Tumor Network - SIGNED_1.1 (tuned)"
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

# Tumor Unsigned Dendrogram (tuned)
plotDendroAndColors(
  brca_netwk_tumor_unsigned_1$dendrograms[[1]],
  colors_tum_unsigned_1[brca_netwk_tumor_unsigned_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Tumor Network - UNSIGNED (tuned)"
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


# normal signed 1 (tuned)
plotDendroAndColors(
  brca_netwk_normal_signed_1_1$dendrograms[[1]],
  colors_norm_signed_1_1[brca_netwk_normal_signed_1_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - SIGNED_1.1 (tuned)"
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

# normal unsigned 1 (tuned)
plotDendroAndColors(
  brca_netwk_normal_unsigned_1_1$dendrograms[[1]],
  colors_norm_unsigned_1_1[brca_netwk_normal_unsigned_1_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Normal Network - UNSIGNED_1.1 (tuned)"
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

# finalised networks:
# brca_netwk_tumor_signed_1, brca_netwk_tumor_unsigned_1
# brca_netwk_normal_signed_1_1, brca_netwk_normal_unsigned_1_1




# MODULE EIGENGENE COMPARISION
# signed and unsigned


# calculating and cleaning eigengens

# Tumor Signed
MEs_tum_signed <- moduleEigengenes(
  t(filtered_tumor), 
  colors = colors_tum_signed_1
)$eigengenes
MEs_tum_signed <- MEs_tum_signed[, colnames(MEs_tum_signed) != "MEgrey"]
MEs_tum_signed <- orderMEs(MEs_tum_signed)
colnames(MEs_tum_signed) <- gsub("ME", "", colnames(MEs_tum_signed))

# Tumor Unsigned
MEs_tum_unsigned <- moduleEigengenes(
  t(filtered_tumor), 
  colors = colors_tum_unsigned_1
)$eigengenes
MEs_tum_unsigned <- MEs_tum_unsigned[, colnames(MEs_tum_unsigned) != "MEgrey"]
MEs_tum_unsigned <- orderMEs(MEs_tum_unsigned)
colnames(MEs_tum_unsigned) <- gsub("ME", "", colnames(MEs_tum_unsigned))

# 3. Normal Signed
MEs_norm_signed <- moduleEigengenes(
  t(filtered_normal), 
  colors = colors_norm_signed_1_1
)$eigengenes
MEs_norm_signed <- MEs_norm_signed[, colnames(MEs_norm_signed) != "MEgrey"]
MEs_norm_signed <- orderMEs(MEs_norm_signed)
colnames(MEs_norm_signed) <- gsub("ME", "", colnames(MEs_norm_signed))

# 4. Normal Unsigned
MEs_norm_unsigned <- moduleEigengenes(
  t(filtered_normal), 
  colors = colors_norm_unsigned_1_1
)$eigengenes
MEs_norm_unsigned <- MEs_norm_unsigned[, colnames(MEs_norm_unsigned) != "MEgrey"]
MEs_norm_unsigned <- orderMEs(MEs_norm_unsigned)
colnames(MEs_norm_unsigned) <- gsub("ME", "", colnames(MEs_norm_unsigned))


# pairwise pearson correlation
cor_tum_signed    <- cor(MEs_tum_signed)
cor_tum_unsigned  <- cor(MEs_tum_unsigned)
cor_norm_signed   <- cor(MEs_norm_signed)
cor_norm_unsigned <- cor(MEs_norm_unsigned)


# verify eigengene dimensions
cat("Tumor Signed Matrix Dimensions \n")
print(dim(cor_tum_signed))

cat("\nTumor Unsigned Matrix Dimensions \n")
print(dim(cor_tum_unsigned))


# EIGENGENE CORRELATION HEATMAPS


library(gplots)

# Step 1: Calculate the 4 internal Pearson correlation matrices
cor_tum_unsigned  <- cor(MEs_tum_unsigned)
cor_tum_signed    <- cor(MEs_tum_signed)
cor_norm_unsigned <- cor(MEs_norm_unsigned)
cor_norm_signed   <- cor(MEs_norm_signed)

# Step 2: Split the plotting window into a 2x2 grid to see all 4 at once
par(mfrow = c(2, 2))

# Tumor Unsigned
heatmap.2(cor_tum_unsigned,
          main = "Tumor Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# Tumor Signed
heatmap.2(cor_tum_signed,
          main = "Tumor Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# Normal Unsigned
heatmap.2(cor_norm_unsigned,
          main = "Normal Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# Normal Signed
heatmap.2(cor_norm_signed,
          main = "Normal Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

# Reset plotting grid back to default single panel
par(mfrow = c(1, 1))


# export to cytoscape format
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
dir.create("modules", showWarnings = FALSE)


# FUNCTION: Reusable pipeline to extract Edge and Node files

export_to_cytoscape <- function(filtered_data, power_val, net_type, colors_vector, prefix) {
  
  cat("Processing TOM for", prefix, "\n")
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
  
  # STRIP DECIMALS FOR DATABASE LOOKUP
  gene1_clean <- gsub("\\..*", "", edge_list$gene1)
  gene2_clean <- gsub("\\..*", "", edge_list$gene2)
  
  edge_list$gene1.name <- mapIds(org.Hs.eg.db, keys = gene1_clean, column = "SYMBOL", keytype = "ENSEMBL")
  edge_list$gene2.name <- mapIds(org.Hs.eg.db, keys = gene2_clean, column = "SYMBOL", keytype = "ENSEMBL")
  
  # Fill in missing symbols with original Ensembl IDs as a fallback safety check
  edge_list$gene1.name <- ifelse(is.na(edge_list$gene1.name), edge_list$gene1, edge_list$gene1.name)
  edge_list$gene2.name <- ifelse(is.na(edge_list$gene2.name), edge_list$gene2, edge_list$gene2.name)
  
  # 5. Construct the corresponding Node List mapping genes to their final colors
  node_list <- data.frame(
    GeneID = row.names(filtered_data),
    ModuleColor = colors_vector
  )
  
  # STRIP DECIMALS FOR NODE LOOKUP
  node_clean <- gsub("\\..*", "", node_list$GeneID)
  
  node_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = node_clean, column = "SYMBOL", keytype = "ENSEMBL")
  node_list$GeneSymbol <- ifelse(is.na(node_list$GeneSymbol), node_list$GeneID, node_list$GeneSymbol)
  
  # 6. Save data tables to your "modules" folder as clean, tab-separated text files
  write.table(edge_list, file = paste0("modules/Edges_", prefix, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(node_list, file = paste0("modules/Nodes_", prefix, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Successfully exported Edge and Node files for:", prefix, "\n\n")
}

# EXECUTE RUNS FOR ALL 4 INDIVIDUAL CONFIGURATIONS


# Run 1: Tumor Unsigned
export_to_cytoscape(filtered_tumor, power_tumor_unsigned, "unsigned", colors_tum_unsigned_1, "Tumor_Unsigned")

# Run 2: Tumor Signed
export_to_cytoscape(filtered_tumor, power_tumor_signed, "signed", colors_tum_signed_1, "Tumor_Signed")

# Run 3: Normal Unsigned
export_to_cytoscape(filtered_normal, power_normal_unsigned_1, "unsigned", colors_norm_unsigned_1_1, "Normal_Unsigned")

# Run 4: Normal Signed
export_to_cytoscape(filtered_normal, power_normal_signed_1, "signed", colors_norm_signed_1_1, "Normal_Signed")

# SEE THE CYTOSCAPE_GUIDE PDF TO UNDERSTAND HOW TO USE THESE FILES TO VISUALISE THE NETWORKS IN CYTOSCAPE


# MODULE PRESERVATION ANALYSIS

library(WGCNA)
library(tidyverse)
allowWGCNAThreads()


# BLOCK 1: SIGNED NETWORK MODULE PRESERVATION

cat("\nSTEP 1: Preparing Inputs for SIGNED Network\n")

# 1. Package the Transposed Data Matrices (Samples as Rows, Genes as Columns)
multiData_signed <- list(
  Tumor  = list(data = t(filtered_tumor)),
  Normal = list(data = t(filtered_normal))
)

# 2. Extract the SIGNED Tumor & Normal Character Color Vectors
tumor_colors_signed  <- labels2colors(brca_netwk_tumor_signed_1$colors)
names(tumor_colors_signed) <- names(brca_netwk_tumor_signed_1$colors)

normal_colors_signed <- labels2colors(brca_netwk_normal_signed_1_1$colors)
names(normal_colors_signed) <- names(brca_netwk_normal_signed_1_1$colors)

multiColor_signed <- list(
  Tumor  = tumor_colors_signed,
  Normal = normal_colors_signed
)

# 3. Quick Integrity Check (Both must return TRUE)
if(all(names(multiColor_signed$Tumor) %in% rownames(filtered_tumor)) && 
   all(names(multiColor_signed$Normal) %in% rownames(filtered_normal))) {
  cat("Sanity checks passed for Signed Network. Proceeding to calculations\n")
} else {
  stop("Error: Gene name mismatch detected in Signed Network structures!")
}

# 4. Clean multiData for Zero Variance / Missing Entries
cat("\nCleaning multiData for zero-variance genes/samples...\n")

# Run the multi-set good samples and genes check
gsg <- goodSamplesGenesMS(multiData_signed, verbose = 3)

if (!gsg$allOK) {
  cat("Removing genes:", sum(!gsg$goodGenes), "\n")
  
  # Manually subset the data using standard R indexing
  # Keeps only the good samples for each set, and the good genes across both
  multiData_signed$Tumor$data  <- multiData_signed$Tumor$data[gsg$goodSamples[[1]], gsg$goodGenes]
  multiData_signed$Normal$data <- multiData_signed$Normal$data[gsg$goodSamples[[2]], gsg$goodGenes]
  
  # Synchronize the color vectors to match the remaining genes
  remaining_genes          <- colnames(multiData_signed$Tumor$data)
  multiColor_signed$Tumor  <- multiColor_signed$Tumor[remaining_genes]
  multiColor_signed$Normal <- multiColor_signed$Normal[remaining_genes]
  
  cat("Cleaned multiData and synchronized color vectors.\n")
} else {
  cat("No zero-variance or missing data issues found.\n")
}

# 5. Run the Permutation Math for Signed
cat("\nRunning SIGNED Module Preservation (100 loops)\n")
preservation_signed <- modulePreservation(
  multiData         = multiData_signed,
  multiColor        = multiColor_signed,
  referenceNetworks = 1,       
  nPermutations     = 100,     
  randomSeed        = 12345,   
  verbose           = 3
)
# 5. Extract Statistics and Filter out the Artificial "Gold" Control Module
stats_signed <- preservation_signed$preservation$Z$ref.Tumor$inColumnsAlsoPresentIn.Normal
stats_signed <- stats_signed[rownames(stats_signed) != "gold", ]

cat("\nSIGNED Network Permutations Complete!\n")
print(head(stats_signed))



# BLOCK 2: UNSIGNED NETWORK MODULE PRESERVATION

cat("\nSTEP 2: Preparing Inputs for UNSIGNED Network\n")

# 1. Package Data Matrices
multiData_unsigned <- list(
  Tumor  = list(data = t(filtered_tumor)),
  Normal = list(data = t(filtered_normal))
)

# 2. Extract the UNSIGNED Tumor & Normal Character Color Vectors
tumor_colors_unsigned  <- labels2colors(brca_netwk_tumor_unsigned_1$colors)
names(tumor_colors_unsigned) <- names(brca_netwk_tumor_unsigned_1$colors)

normal_colors_unsigned <- labels2colors(brca_netwk_normal_unsigned_1_1$colors)
names(normal_colors_unsigned) <- names(brca_netwk_normal_unsigned_1_1$colors)

multiColor_unsigned <- list(
  Tumor  = tumor_colors_unsigned,
  Normal = normal_colors_unsigned
)

# 3. Quick Integrity Check (Both must return TRUE)
if(all(names(multiColor_unsigned$Tumor) %in% rownames(filtered_tumor)) && 
   all(names(multiColor_unsigned$Normal) %in% rownames(filtered_normal))) {
  cat("Sanity checks passed for Unsigned Network. Proceeding to calculations...\n")
} else {
  stop("Error: Gene name mismatch detected in Unsigned Network structures!")
}

# 4. Clean multiData for Zero Variance / Missing Entries
cat("\nCleaning multiData for zero-variance genes/samples\n")

# Run the multi-set good samples and genes check
gsg <- goodSamplesGenesMS(multiData_unsigned, verbose = 3)

if (!gsg$allOK) {
  cat("Removing genes:", sum(!gsg$goodGenes), "\n")
  
  # Manually subset the data using standard R indexing
  # Keeps only the good samples for each set, and the good genes across both
  multiData_unsigned$Tumor$data  <- multiData_unsigned$Tumor$data[gsg$goodSamples[[1]], gsg$goodGenes]
  multiData_unsigned$Normal$data <- multiData_unsigned$Normal$data[gsg$goodSamples[[2]], gsg$goodGenes]
  
  # Synchronize the color vectors to match the remaining genes
  remaining_genes          <- colnames(multiData_unsigned$Tumor$data)
  multiColor_unsigned$Tumor  <- multiColor_unsigned$Tumor[remaining_genes]
  multiColor_unsigned$Normal <- multiColor_unsigned$Normal[remaining_genes]
  
  cat("Cleaned multiData and synchronized color vectors.\n")
} else {
  cat("No zero-variance or missing data issues found.\n")
}

# 5. Run the Permutation Math for Signed
cat("\nRunning UNSIGNED Module Preservation (100 loops)\n")
preservation_unsigned <- modulePreservation(
  multiData         = multiData_unsigned,
  multiColor        = multiColor_unsigned,
  referenceNetworks = 1,       
  nPermutations     = 100,     
  randomSeed        = 12345,   
  verbose           = 3
)
# 5. Extract Statistics and Filter out the Artificial "Gold" Control Module
stats_unsigned <- preservation_unsigned$preservation$Z$ref.Tumor$inColumnsAlsoPresentIn.Normal
stats_unsigned <- stats_unsigned[rownames(stats_unsigned) != "gold", ]

cat("\nUNSIGNED Network Permutations Complete!\n")
print(head(stats_unsigned))


# EXTRACT STATS & GENERATE Z-SUMMARY DOT PLOTS

library(ggplot2)


# PART A: Plotting the SIGNED Network Preservation

# Create a clean dataframe for ggplot
plot_data_signed <- data.frame(
  Module    = rownames(stats_signed),
  Z_summary = as.numeric(stats_signed$Zsummary.pres)
)

# Generate the Dot Plot for Signed
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


# PART B: Plotting the UNSIGNED Network Preservation

# Create a clean dataframe for ggplot
plot_data_unsigned <- data.frame(
  Module    = rownames(stats_unsigned),
  Z_summary = as.numeric(stats_unsigned$Zsummary.pres)
)

# Generate the Dot Plot for Unsigned
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



# GO ENRICHMENT

library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)


# Clean out stale empty CSVs from the failed previous run

if (dir.exists("enrich/GO_Signed"))   file.remove(list.files("enrich/GO_Signed",   full.names = TRUE))
if (dir.exists("enrich/GO_Unsigned")) file.remove(list.files("enrich/GO_Unsigned", full.names = TRUE))

dir.create("enrich",             showWarnings = FALSE)
dir.create("enrich/GO_Signed",   showWarnings = FALSE)
dir.create("enrich/GO_Unsigned", showWarnings = FALSE)


# STEP 1: Core GO Enrichment Function
# strips Ensembl version decimals (e.g. ENSG00000000005.6 -> ENSG00000000005)
# before passing to org.Hs.eg.db, which only accepts bare Ensembl IDs

perform_go_enrichment <- function(gene_list, ontology, output_path) {
  
  # Strip version suffixes
  gene_list_clean <- gsub("\\..*", "", gene_list)
  
  go_results <- enrichGO(
    gene          = gene_list_clean,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Only write CSV if results exist (avoids empty file crash downstream)
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    write.csv(as.data.frame(go_results), file = output_path, row.names = TRUE)
  }
  
  return(go_results)
}


# STEP 2: Plot generation function for one module 

plot_go_for_module <- function(go_bp, go_cc, go_mf, mod_name, network_type, out_dir) {
  
  plots <- list()
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    plots[["BP"]] <- dotplot(go_bp, showCategory = 10, font.size = 9, label_format = 60) +
      scale_size_continuous(range = c(2, 7)) +
      theme_minimal() +
      ggtitle(paste0("GO Biological Process (", mod_name, " - ", network_type, ")"))
  }
  
  if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
    plots[["CC"]] <- dotplot(go_cc, showCategory = 10, font.size = 9, label_format = 60) +
      scale_size_continuous(range = c(2, 7)) +
      theme_minimal() +
      ggtitle(paste0("GO Cellular Component (", mod_name, " - ", network_type, ")"))
  }
  
  if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
    plots[["MF"]] <- dotplot(go_mf, showCategory = 10, font.size = 9, label_format = 60) +
      scale_size_continuous(range = c(2, 7)) +
      theme_minimal() +
      ggtitle(paste0("GO Molecular Function (", mod_name, " - ", network_type, ")"))
  }
  
  if (length(plots) > 0) {
    combined <- ggarrange(plotlist = plots, ncol = 1, nrow = length(plots))
    out_file <- paste0(out_dir, "/GO_combined_", mod_name, ".png")
    ggsave(out_file, combined, width = 10, height = 5 * length(plots), dpi = 300)
    cat("  Saved combined dotplot for", mod_name, "to", out_file, "\n")
  } else {
    cat("  No significant GO terms found for any ontology in module:", mod_name, "\n")
  }
}


# STEP 3: SIGNED Network GO Enrichment

cat("\nSIGNED NETWORK: Building Gene Lists\n")


module_names_signed   <- unique(tumor_colors_signed)
genes_in_modules_signed <- lapply(module_names_signed, function(mod) {
  names(tumor_colors_signed[tumor_colors_signed == mod])
})
names(genes_in_modules_signed) <- module_names_signed

signed_modules_to_test <- module_names_signed[module_names_signed != "grey"]
cat("Modules to test (Signed):", paste(signed_modules_to_test, collapse = ", "), "\n")

cat("\nRunning GO Enrichment for Signed Modules\n")

go_results_signed <- list()

for (mod in signed_modules_to_test) {
  cat("\nProcessing module:", mod, "(Signed)\n")
  gene_list <- genes_in_modules_signed[[mod]]
  go_results_signed[[mod]] <- list()
  
  for (ont in c("BP", "CC", "MF")) {
    out_path <- paste0("enrich/GO_Signed/GO_", ont, "_", mod, ".csv")
    cat("  Running", ont, "...\n")
    
    result <- tryCatch(
      perform_go_enrichment(gene_list, ont, out_path),
      error = function(e) {
        cat("  No significant terms for Signed:", mod, "(", ont, ") -", e$message, "\n")
        return(NULL)
      }
    )
    go_results_signed[[mod]][[ont]] <- result
  }
  
  # Generate combined 3-panel dotplot for this module
  plot_go_for_module(
    go_bp       = go_results_signed[[mod]][["BP"]],
    go_cc       = go_results_signed[[mod]][["CC"]],
    go_mf       = go_results_signed[[mod]][["MF"]],
    mod_name    = mod,
    network_type = "Signed",
    out_dir     = "enrich/GO_Signed"
  )
}

cat("\nSIGNED GO Enrichment Complete\n")

# FOUND NO GO ENRICHMENT RESULTS FOR THE RED MODULE

# Counting how many clean genes we actually have in the red module
red_genes <- names(tumor_colors_signed[tumor_colors_signed == "red"])
red_clean <- gsub("\\..*", "", red_genes)
print(length(red_clean))

# Checking how many of them actually exist in the database
library(org.Hs.eg.db)
mapped <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys    = red_clean, 
  columns = "SYMBOL", 
  keytype = "ENSEMBL"
)

# Seeing how many actually mapped successfully
successful_maps <- na.omit(mapped)
cat("Total genes in Red module:", length(red_clean), "\n")
cat("Successfully mapped to Symbols:", nrow(successful_maps), "\n")

# STEP 4: UNSIGNED Network GO Enrichment

cat("\nUNSIGNED NETWORK: Building Gene Lists\n")

module_names_unsigned   <- unique(tumor_colors_unsigned)
genes_in_modules_unsigned <- lapply(module_names_unsigned, function(mod) {
  names(tumor_colors_unsigned[tumor_colors_unsigned == mod])
})
names(genes_in_modules_unsigned) <- module_names_unsigned

unsigned_modules_to_test <- module_names_unsigned[module_names_unsigned != "grey"]
cat("Modules to test (Unsigned):", paste(unsigned_modules_to_test, collapse = ", "), "\n")

cat("\nRunning GO Enrichment for Unsigned Modules\n")

go_results_unsigned <- list()

for (mod in unsigned_modules_to_test) {
  cat("\nProcessing module:", mod, "(Unsigned)\n")
  gene_list <- genes_in_modules_unsigned[[mod]]
  go_results_unsigned[[mod]] <- list()
  
  for (ont in c("BP", "CC", "MF")) {
    out_path <- paste0("enrich/GO_Unsigned/GO_", ont, "_", mod, ".csv")
    cat("  Running", ont, "...\n")
    
    result <- tryCatch(
      perform_go_enrichment(gene_list, ont, out_path),
      error = function(e) {
        cat("  No significant terms for Unsigned:", mod, "(", ont, ") -", e$message, "\n")
        return(NULL)
      }
    )
    go_results_unsigned[[mod]][[ont]] <- result
  }
  
  # Generate combined 3-panel dotplot for this module
  plot_go_for_module(
    go_bp        = go_results_unsigned[[mod]][["BP"]],
    go_cc        = go_results_unsigned[[mod]][["CC"]],
    go_mf        = go_results_unsigned[[mod]][["MF"]],
    mod_name     = mod,
    network_type = "Unsigned",
    out_dir      = "enrich/GO_Unsigned"
  )
}

cat("\nUNSIGNED GO Enrichment Complete\n")


# STEP 5: Summary

cat("\nSUMMARY\n")

cat("\nSIGNED modules with at least one significant GO term:\n")
for (mod in signed_modules_to_test) {
  has_results <- sapply(c("BP","CC","MF"), function(ont) {
    r <- go_results_signed[[mod]][[ont]]
    !is.null(r) && nrow(as.data.frame(r)) > 0
  })
  cat(" ", mod, "->", paste(names(has_results)[has_results], collapse = ", "),
      if (any(has_results)) "" else "(no significant terms)", "\n")
}

cat("\nUNSIGNED modules with at least one significant GO term:\n")
for (mod in unsigned_modules_to_test) {
  has_results <- sapply(c("BP","CC","MF"), function(ont) {
    r <- go_results_unsigned[[mod]][[ont]]
    !is.null(r) && nrow(as.data.frame(r)) > 0
  })
  cat(" ", mod, "->", paste(names(has_results)[has_results], collapse = ", "),
      if (any(has_results)) "" else "(no significant terms)", "\n")
}

# HUB GENE IDENTIFICATION

# Hub genes are the most highly connected members of each co-expression module.
# They are identified by computing Module Membership (kME): the Pearson
# correlation between each gene's expression profile and the module eigengene.
# Genes with the highest absolute kME are the hub genes.


library(WGCNA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Number of top hub genes to extract per module
N_HUBS <- 20

dir.create("hub_genes", showWarnings = FALSE)

# SECTION 1: SIGNED NETWORK HUB GENES


cat("\nSIGNED NETWORK: Computing Module Membership\n")


MM_signed <- as.data.frame(
  cor(t(filtered_tumor), MEs_tum_signed, use = "pairwise.complete.obs")
)

# The column names of MM_signed now match the module colour names in MEs_tum_signed.
# Confirm:
cat("Module membership matrix dimensions:", dim(MM_signed), "\n")
cat("Modules covered:", colnames(MM_signed), "\n")


# Extract top N hub genes per module (signed)


# Get the module colour names (excluding grey, which was already removed from MEs)
signed_modules <- colnames(MM_signed)

hub_list_signed <- list()

for (mod in signed_modules) {
  
  # kME scores for this module across all genes
  kme_scores <- MM_signed[[mod]]
  names(kme_scores) <- rownames(filtered_tumor)
  
  # Sort by absolute kME descending and take top N
  top_genes_ensembl <- names(sort(abs(kme_scores), decreasing = TRUE))[1:N_HUBS]
  top_kme_values    <- kme_scores[top_genes_ensembl]
  
  # Strip Ensembl version suffixes for database lookup
  top_genes_clean <- gsub("\\..*", "", top_genes_ensembl)
  
  # Translate to HGNC gene symbols
  gene_symbols <- tryCatch(
    mapIds(org.Hs.eg.db,
           keys      = top_genes_clean,
           column    = "SYMBOL",
           keytype   = "ENSEMBL",
           multiVals = "first"),
    error = function(e) setNames(top_genes_ensembl, top_genes_ensembl)
  )
  
  hub_list_signed[[mod]] <- data.frame(
    Module        = mod,
    EnsemblID     = top_genes_ensembl,
    GeneSymbol    = gene_symbols,
    kME           = round(top_kme_values, 4),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
  
  cat("  Module", mod, "-- top hub gene:", gene_symbols[1],
      "(kME =", round(top_kme_values[1], 3), ")\n")
}

# Combine into a single data frame
hub_genes_signed <- do.call(rbind, hub_list_signed)
rownames(hub_genes_signed) <- NULL

# Save to CSV
write.csv(hub_genes_signed,
          "hub_genes/hub_genes_signed_top20.csv",
          row.names = FALSE)

cat("\nSigned hub gene table saved to hub_genes/hub_genes_signed_top20.csv\n")
print(hub_genes_signed)



# SECTION 2: UNSIGNED NETWORK HUB GENES


cat("\nUNSIGNED NETWORK: Computing Module Membership\n")

MM_unsigned <- as.data.frame(
  cor(t(filtered_tumor), MEs_tum_unsigned, use = "pairwise.complete.obs")
)

cat("Module membership matrix dimensions:", dim(MM_unsigned), "\n")
cat("Modules covered:", colnames(MM_unsigned), "\n")

unsigned_modules <- colnames(MM_unsigned)

hub_list_unsigned <- list()

for (mod in unsigned_modules) {
  
  kme_scores <- MM_unsigned[[mod]]
  names(kme_scores) <- rownames(filtered_tumor)
  
  top_genes_ensembl <- names(sort(abs(kme_scores), decreasing = TRUE))[1:N_HUBS]
  top_kme_values    <- kme_scores[top_genes_ensembl]
  
  top_genes_clean <- gsub("\\..*", "", top_genes_ensembl)
  
  gene_symbols <- tryCatch(
    mapIds(org.Hs.eg.db,
           keys      = top_genes_clean,
           column    = "SYMBOL",
           keytype   = "ENSEMBL",
           multiVals = "first"),
    error = function(e) setNames(top_genes_ensembl, top_genes_ensembl)
  )
  
  hub_list_unsigned[[mod]] <- data.frame(
    Module        = mod,
    EnsemblID     = top_genes_ensembl,
    GeneSymbol    = gene_symbols,
    kME           = round(top_kme_values, 4),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
  
  cat("  Module", mod, "-- top hub gene:", gene_symbols[1],
      "(kME =", round(top_kme_values[1], 3), ")\n")
}

hub_genes_unsigned <- do.call(rbind, hub_list_unsigned)
rownames(hub_genes_unsigned) <- NULL

write.csv(hub_genes_unsigned,
          "hub_genes/hub_genes_unsigned_top20.csv",
          row.names = FALSE)

cat("\nUnsigned hub gene table saved to hub_genes/hub_genes_unsigned_top20.csv\n")
print(hub_genes_unsigned)



# SECTION 3: VISUALISE kME DISTRIBUTIONS (bar charts per module)
# Shows the kME score of the top N hub genes for each module in both networks.


cat("\nGenerating kME bar plots\n")

plot_hub_bar <- function(hub_df, network_type) {
  
  modules_to_plot <- unique(hub_df$Module)
  
  for (mod in modules_to_plot) {
    
    df_mod <- hub_df[hub_df$Module == mod, ]
    
    # Use GeneSymbol if available, otherwise fall back to EnsemblID
    df_mod$Label <- ifelse(
      is.na(df_mod$GeneSymbol) | df_mod$GeneSymbol == "",
      df_mod$EnsemblID,
      df_mod$GeneSymbol
    )
    
    # Sort by absolute kME for plotting
    df_mod <- df_mod[order(abs(df_mod$kME), decreasing = FALSE), ]
    df_mod$Label <- factor(df_mod$Label, levels = df_mod$Label)
    
    p <- ggplot(df_mod, aes(x = abs(kME), y = Label)) +
      geom_bar(stat = "identity", fill = mod, colour = "grey30",
               linewidth = 0.3, alpha = 0.85) +
      labs(
        title    = paste0("Top ", N_HUBS, " Hub Genes -- ",
                          toupper(mod), " Module (", network_type, ")"),
        subtitle = paste0("Ranked by |kME| (Module Membership)"),
        x        = "Module Membership |kME|",
        y        = "Gene Symbol"
      ) +
      xlim(0, 1) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y  = element_text(size = 9),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40")
      )
    
    out_file <- paste0("hub_genes/hub_bar_", network_type, "_", mod, ".png")
    ggsave(out_file, p, width = 7, height = 6, dpi = 300)
    cat("  Saved:", out_file, "\n")
  }
}

plot_hub_bar(hub_genes_signed,   "signed")
plot_hub_bar(hub_genes_unsigned, "unsigned")

cat("\nAll hub gene bar plots saved to hub_genes/\n")



# SECTION 4: COMBINED kME SUMMARY TABLE
# A single flat table covering both networks for easy comparison.

hub_genes_signed$Network   <- "Signed"
hub_genes_unsigned$Network <- "Unsigned"

hub_genes_all <- rbind(hub_genes_signed, hub_genes_unsigned)
hub_genes_all <- hub_genes_all[, c("Network", "Module", "GeneSymbol",
                                   "EnsemblID", "kME")]

write.csv(hub_genes_all,
          "hub_genes/hub_genes_combined_both_networks.csv",
          row.names = FALSE)

cat("\nCombined hub gene table (both networks) saved to:\n")
cat("  hub_genes/hub_genes_combined_both_networks.csv\n")


# summary


cat("HUB GENE IDENTIFICATION COMPLETE\n")


cat("\nSIGNED NETWORK -- Top hub gene per module:\n")
for (mod in unique(hub_genes_signed$Module)) {
  top <- hub_genes_signed[hub_genes_signed$Module == mod, ][1, ]
  cat(sprintf("  %-12s : %-10s (kME = %.3f)\n",
              mod, top$GeneSymbol, top$kME))
}

cat("\nUNSIGNED NETWORK -- Top hub gene per module:\n")
for (mod in unique(hub_genes_unsigned$Module)) {
  top <- hub_genes_unsigned[hub_genes_unsigned$Module == mod, ][1, ]
  cat(sprintf("  %-12s : %-10s (kME = %.3f)\n",
              mod, top$GeneSymbol, top$kME))
}

cat("\nOutput files:\n")
cat("  hub_genes/hub_genes_signed_top20.csv\n")
cat("  hub_genes/hub_genes_unsigned_top20.csv\n")
cat("  hub_genes/hub_genes_combined_both_networks.csv\n")
cat("  hub_genes/hub_bar_signed_<module>.png  (one per module)\n")
cat("  hub_genes/hub_bar_unsigned_<module>.png (one per module)\n")

