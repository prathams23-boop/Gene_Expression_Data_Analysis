# ==============================================================================
#
#         WGCNA PIPELINE -- UNIVERSAL TEMPLATE SCRIPT
#
#  This script is a complete, reusable template for Weighted Gene Co-expression
#  Network Analysis (WGCNA) on bulk RNA-seq data. It covers the full workflow:
#
#    1. Data loading and exploration
#    2. Metadata construction and sanity checks
#    3. DESeq2 normalisation (design = ~1, intercept-only)
#    4. Variance-stabilising transformation (VST)
#    5. Batch effect check via PCA
#    6. Variance-based gene filtering (top 5%) and gene-set alignment
#    7. Soft-thresholding power selection (signed and unsigned)
#    8. Scale independence and mean connectivity diagnostic plots
#    9. Co-expression network construction (signed and unsigned)
#   10. Cluster dendrogram visualisation
#   11. Module eigengene computation and correlation heatmaps
#   12. Network export to Cytoscape (edge lists and node lists)
#   13. Module preservation analysis (permutation test, Z-summary)
#   14. GO enrichment analysis for all modules (BP, CC, MF)
#   15. Individual module GO dotplots
#
#  HOW TO ADAPT THIS TEMPLATE:
#  Search for every line that begins with "# EDIT:" and update it for your
#  specific dataset. Everything else can usually be left unchanged.
#
#  REQUIRED R PACKAGES (install once via Bioconductor and CRAN):
#    install.packages(c("WGCNA","tidyverse","ggplot2","ggpubr","gplots",
#                       "dendextend","VennDiagram"))
#    BiocManager::install(c("DESeq2","genefilter","clusterProfiler",
#                           "org.Hs.eg.db","GO.db","sva"))
#    # Replace org.Hs.eg.db with the correct organism annotation package:
#    #   Human  : org.Hs.eg.db
#    #   Mouse  : org.Mm.eg.db
#    #   Rat    : org.Rn.eg.db
#    #   Zebrafish: org.Dr.eg.db
#
# ==============================================================================


# ==============================================================================
# SECTION 0: USER-EDITABLE CONFIGURATION
# All dataset-specific settings are collected here so you only need to edit
# this block when working with a new dataset.
# ==============================================================================

# ---- File paths --------------------------------------------------------------
# EDIT: Set the full paths to your two CSV gene expression matrices.
# The CSV files must have genes as rows (with Ensembl IDs as row names)
# and samples as columns.

PATH_REFERENCE <- "path/to/your/reference_condition_data.csv"   # e.g. Tumor
PATH_TEST      <- "path/to/your/test_condition_data.csv"         # e.g. Normal

# ---- Condition labels --------------------------------------------------------
# EDIT: Change these to match your experimental conditions.
LABEL_REFERENCE <- "Tumor"    # Label used in metadata for the reference group
LABEL_TEST      <- "Normal"   # Label used in metadata for the test/baseline group

# ---- Sample counts -----------------------------------------------------------
# EDIT: Set the exact number of samples in each matrix.
N_REFERENCE <- 337   # Number of columns in reference data (e.g. 337 tumour samples)
N_TEST      <- 32    # Number of columns in test data (e.g. 32 normal samples)

# ---- Organism annotation database --------------------------------------------
# EDIT: Replace with the correct annotation package for your organism.
# See the header comment for common choices.
ORG_DB_PACKAGE <- "org.Hs.eg.db"

# ---- Soft-thresholding powers ------------------------------------------------
# These are set AFTER running Section 5 (pickSoftThreshold).
# Leave as NA now; fill them in once you have seen the diagnostic plots.
# EDIT: Enter your validated powers here after running Section 5.
POWER_REF_SIGNED    <- NA   # e.g. 8
POWER_REF_UNSIGNED  <- NA   # e.g. 4
POWER_TEST_SIGNED   <- NA   # e.g. 20
POWER_TEST_UNSIGNED <- NA   # e.g. 14

# ---- Network construction parameters -----------------------------------------
# These defaults are appropriate for most RNA-seq datasets.
# EDIT: Adjust deepSplit (0-4) or minModuleSize if your modules are too large
# or too fragmented.
DEEP_SPLIT       <- 2     # Sensitivity of tree-cutting (0=coarse, 4=fine)
MIN_MODULE_SIZE  <- 30    # Minimum number of genes in a valid module
MERGE_CUT_HEIGHT <- 0.25  # Merge modules with eigengene correlation > 1-0.25 = 0.75

# ---- Module preservation parameters -----------------------------------------
N_PERMUTATIONS <- 100     # 100 is standard; use 200+ for publication
RANDOM_SEED    <- 12345   # Fix for reproducibility

# ---- Output directories ------------------------------------------------------
# EDIT: Change if you want outputs in a different location.
DIR_MODULES <- "modules"       # Cytoscape edge/node files
DIR_ENRICH  <- "enrich/GO"     # GO enrichment CSVs and plots
DIR_FIGURES <- "figures"       # Saved PDF/PNG figures

# ---- Cytoscape edge weight threshold -----------------------------------------
# Only gene pairs with TOM > this value are written to the edge list.
# Lower = more edges (denser network); higher = fewer edges (sparser).
TOM_THRESHOLD <- 0.1

# ---- GO enrichment parameters ------------------------------------------------
GO_PVAL_CUTOFF  <- 0.05
GO_QVAL_CUTOFF  <- 0.05
GO_TOP_CATEGORY <- 10     # Number of top terms shown in each dotplot


# ==============================================================================
# SECTION 1: LOAD LIBRARIES
# ==============================================================================

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gplots)
library(dendextend)
library(VennDiagram)
library(dplyr)
library(DESeq2)
library(genefilter)
library(clusterProfiler)
library(GO.db)

# Load the organism annotation database defined in Section 0
library(ORG_DB_PACKAGE, character.only = TRUE)
org_db <- get(ORG_DB_PACKAGE)

# Create output directories
dir.create(DIR_MODULES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_ENRICH,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

cat("Libraries loaded and output directories created.\n")


# ==============================================================================
# SECTION 2: DATA LOADING AND INITIAL EXPLORATION
# ==============================================================================

cat("Loading data...\n")

# row.names=1 tells R that the first column of the CSV contains gene identifiers
# (typically Ensembl IDs) and should be used as row names, not treated as data.
reference_data <- read.csv(PATH_REFERENCE, row.names = 1)
test_data      <- read.csv(PATH_TEST,      row.names = 1)

# ---- Inspect dimensions -------------------------------------------------------
cat("Reference data dimensions:", dim(reference_data), "\n")
cat("Test data dimensions:     ", dim(test_data), "\n")

# ---- Quick preview of first 5 rows and columns --------------------------------
cat("\nFirst 5 rows and columns of reference data:\n")
print(reference_data[1:5, 1:5])

cat("\nFirst 5 rows and columns of test data:\n")
print(test_data[1:5, 1:5])

# ---- Optional: check gene ID format -------------------------------------------
cat("\nSample gene IDs from reference data:\n")
print(head(rownames(reference_data)))


# ==============================================================================
# SECTION 3: METADATA CONSTRUCTION AND SANITY CHECKS
# ==============================================================================

# DESeq2 requires a metadata data frame whose row names exactly match the
# column names of the count matrix. The Condition column labels each sample.

metadata_reference <- data.frame(
  Sample    = colnames(reference_data),
  Condition = rep(LABEL_REFERENCE, N_REFERENCE)
)
rownames(metadata_reference) <- metadata_reference$Sample

metadata_test <- data.frame(
  Sample    = colnames(test_data),
  Condition = rep(LABEL_TEST, N_TEST)
)
rownames(metadata_test) <- metadata_test$Sample

# ---- Sanity checks: both must return TRUE ------------------------------------
check_ref  <- all(colnames(reference_data) == metadata_reference$Sample)
check_test <- all(colnames(test_data) == metadata_test$Sample)

cat("Reference metadata check:", check_ref,  "\n")
cat("Test metadata check:     ", check_test, "\n")

if (!check_ref || !check_test) {
  stop("FATAL: Sample name mismatch between data matrix and metadata. ",
       "Check that your CSV column names match the metadata exactly.")
}


# ==============================================================================
# SECTION 4: DESEQ2 NORMALISATION (design = ~1)
# ==============================================================================
# The intercept-only design (~1) tells DESeq2 to estimate size factors for
# depth correction WITHOUT fitting a differential expression model. This is
# essential for WGCNA: we want depth-corrected counts that are NOT pre-adjusted
# for condition differences, so the network reflects true co-expression patterns
# rather than an imposed group contrast.
#
# DESeq2 will warn: "the design is ~1 (just an intercept). Is this intended?"
# This warning is expected and safe to ignore.

cat("Running DESeq2 normalisation for reference data...\n")

dds_reference <- DESeqDataSetFromMatrix(
  countData = round(reference_data),  # round() ensures integer counts
  colData   = metadata_reference,
  design    = ~ 1                     # Intercept only; no group comparison
)
dds_reference <- DESeq(dds_reference)
normalized_reference <- counts(dds_reference, normalized = TRUE)

cat("Running DESeq2 normalisation for test data...\n")

dds_test <- DESeqDataSetFromMatrix(
  countData = round(test_data),
  colData   = metadata_test,
  design    = ~ 1
)
dds_test <- DESeq(dds_test)
normalized_test <- counts(dds_test, normalized = TRUE)

cat("Normalisation complete.\n")


# ==============================================================================
# SECTION 5: VARIANCE-STABILISING TRANSFORMATION (VST)
# ==============================================================================
# VST corrects the mean-variance dependence inherent in RNA-seq count data.
# Without VST, highly expressed genes dominate the correlation landscape simply
# because their larger counts produce larger absolute variance. VST makes the
# variance approximately constant (homoskedastic) across all expression levels,
# giving every gene an equal, unbiased opportunity to contribute to the
# co-expression network.

cat("Applying Variance-Stabilising Transformation...\n")

vsd_reference <- varianceStabilizingTransformation(dds_reference)
vsd_test      <- varianceStabilizingTransformation(dds_test)


# ==============================================================================
# SECTION 6: BATCH EFFECT CHECK (PCA)
# ==============================================================================
# TCGA data is collected across multiple tissue source sites (TSS) with
# different equipment and protocols. A PCA coloured by batch proxy (e.g. TSS
# code from the TCGA barcode) reveals whether technical artefacts dominate the
# variance structure.
#
# TCGA barcode format: TCGA-XX-XXXX-XXX  (positions 6-7 are the TSS code)
# Adjust substr() indices if your barcode format is different.
#
# IF YOU SEE DISTINCT COLOUR CLUSTERS in the PCA (Scenario B), you MUST run
# ComBat correction (from the sva package) before proceeding to gene filtering.
# If the dots look like random mixed confetti (Scenario A), proceed directly.

cat("Generating PCA batch-effect check plot for reference data...\n")

# EDIT: Adjust the substr indices to extract a batch proxy from your barcodes.
# For TCGA barcodes (TCGA.XX.XXXX format), characters 6-7 give the TSS code.
tss_reference <- substr(colnames(vsd_reference), 6, 7)

pca_reference <- prcomp(t(assay(vsd_reference)))
pca_df_ref    <- data.frame(
  PC1   = pca_reference$x[, 1],
  PC2   = pca_reference$x[, 2],
  Batch = tss_reference
)

p_pca_ref <- ggplot(pca_df_ref, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2.5, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = paste0("Batch Effect Check: ", LABEL_REFERENCE, " PCA"),
    x = paste0("PC1: ",
               round(summary(pca_reference)$importance[2, 1] * 100, 1), "% variance"),
    y = paste0("PC2: ",
               round(summary(pca_reference)$importance[2, 2] * 100, 1), "% variance")
  )

ggsave(file.path(DIR_FIGURES, "pca_batch_reference.pdf"),
       p_pca_ref, width = 7, height = 5)

print(p_pca_ref)
cat("PCA plot saved. Inspect it before proceeding.\n")
cat("Scenario A (safe): dots form a random mixed cloud.\n")
cat("Scenario B (action needed): dots cluster by colour = batch effect detected.\n")
cat("If Scenario B: run ComBat from the sva package before continuing.\n\n")

# ---- Optional ComBat correction (uncomment if batch effect is detected) ------
# library(sva)
# batch_vector <- as.numeric(as.factor(tss_reference))
# assay(vsd_reference) <- ComBat(dat    = assay(vsd_reference),
#                                 batch  = batch_vector,
#                                 mod    = NULL,
#                                 par.prior = TRUE)
# cat("ComBat correction applied to reference data.\n")


# ==============================================================================
# SECTION 7: VARIANCE-BASED GENE FILTERING AND ALIGNMENT
# ==============================================================================
# After VST, keep only the top 5% most variable genes in each dataset.
# This removes flatline housekeeping genes and low-expression noise that
# carry no co-expression information. The filter is variance-based (NOT
# differential-expression-based): pre-selecting DE genes would artificially
# inflate module-condition correlations.
#
# After independent filtering, the two gene sets must be aligned to their
# intersection so that module preservation analysis can compare networks
# built on identical gene universes.

cat("Computing per-gene variance and filtering to top 5%...\n")

rv_reference <- rowVars(assay(vsd_reference))
rv_test      <- rowVars(assay(vsd_test))

q95_reference <- quantile(rv_reference, 0.95)
q95_test      <- quantile(rv_test,      0.95)

filtered_reference <- assay(vsd_reference)[rv_reference > q95_reference, ]
filtered_test      <- assay(vsd_test)[rv_test > q95_test, ]

cat("Genes retained (reference, before alignment):", nrow(filtered_reference), "\n")
cat("Genes retained (test, before alignment):     ", nrow(filtered_test), "\n")

# ---- Align to common gene set ------------------------------------------------
common_genes       <- intersect(rownames(filtered_reference), rownames(filtered_test))
filtered_reference <- filtered_reference[common_genes, ]
filtered_test      <- filtered_test[common_genes, ]

cat("Genes after alignment (common to both):", length(common_genes), "\n")
cat("Final reference matrix:", dim(filtered_reference), "\n")
cat("Final test matrix:     ", dim(filtered_test), "\n")

if (!identical(dim(filtered_reference)[1], dim(filtered_test)[1])) {
  stop("Gene alignment failed: row counts still differ. Check rownames.")
}


# ==============================================================================
# SECTION 8: SOFT-THRESHOLDING POWER SELECTION
# ==============================================================================
# WGCNA raises the Pearson correlation between each gene pair to a power beta
# to create an adjacency matrix that approximates scale-free topology.
# The goal is to find the LOWEST power at which the scale-free topology fit
# R^2 >= 0.90 AND the slope of the connectivity distribution is negative.
#
# This section runs pickSoftThreshold() for all four configurations.
# AFTER examining the output tables, enter your chosen powers in Section 0
# (POWER_REF_SIGNED, POWER_REF_UNSIGNED, POWER_TEST_SIGNED, POWER_TEST_UNSIGNED).

allowWGCNAThreads()  # Use all available CPU cores

powers <- c(1:10, seq(from = 10, to = 20, by = 2))

cat("Running pickSoftThreshold for reference SIGNED network...\n")
sft_ref_signed <- pickSoftThreshold(
  t(filtered_reference),
  powerVector = powers,
  networkType = "signed",
  verbose     = 3
)

cat("Running pickSoftThreshold for reference UNSIGNED network...\n")
sft_ref_unsigned <- pickSoftThreshold(
  t(filtered_reference),
  powerVector = powers,
  networkType = "unsigned",
  verbose     = 3
)

cat("Running pickSoftThreshold for test SIGNED network...\n")
sft_test_signed <- pickSoftThreshold(
  t(filtered_test),
  powerVector = powers,
  networkType = "signed",
  verbose     = 3
)

cat("Running pickSoftThreshold for test UNSIGNED network...\n")
sft_test_unsigned <- pickSoftThreshold(
  t(filtered_test),
  powerVector = powers,
  networkType = "unsigned",
  verbose     = 3
)

# ---- Print power estimates ---------------------------------------------------
cat("\n--- Soft-thresholding power estimates ---\n")
cat("Reference Signed:    ", sft_ref_signed$powerEstimate,   "\n")
cat("Reference Unsigned:  ", sft_ref_unsigned$powerEstimate, "\n")
cat("Test Signed:         ", sft_test_signed$powerEstimate,  "\n")
cat("Test Unsigned:       ", sft_test_unsigned$powerEstimate,"\n")

cat("\nNOW: Examine the tables above and the plots below.\n")
cat("Enter your chosen powers in Section 0 (POWER_* variables), then re-run.\n\n")


# ==============================================================================
# SECTION 9: SCALE INDEPENDENCE AND MEAN CONNECTIVITY DIAGNOSTIC PLOTS
# ==============================================================================
# These plots visualise the pickSoftThreshold() output. They are required by
# reviewers to justify the chosen soft-thresholding powers.
#
# Scale independence plot: red numbers should form an upside-down L-shape.
#   The first number at or above the horizontal red line (R^2 = 0.90) is
#   your winning power.
# Mean connectivity plot:  numbers should follow a smooth exponential decay.
#   The selected power should sit at the 'knee' of this decay.

cex1 <- 0.85  # Text size for power number labels

# ---- Reference plots ---------------------------------------------------------
pdf(file.path(DIR_FIGURES, "sft_reference.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

# Scale independence -- reference signed
plot(sft_ref_signed$fitIndices[, 1],
     -sign(sft_ref_signed$fitIndices[, 3]) * sft_ref_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = paste0("Scale independence (", LABEL_REFERENCE, " Signed)"))
text(sft_ref_signed$fitIndices[, 1],
     -sign(sft_ref_signed$fitIndices[, 3]) * sft_ref_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

# Mean connectivity -- reference signed
plot(sft_ref_signed$fitIndices[, 1],
     sft_ref_signed$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste0("Mean connectivity (", LABEL_REFERENCE, " Signed)"))
text(sft_ref_signed$fitIndices[, 1],
     sft_ref_signed$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

dev.off()

# ---- Reference unsigned plots ------------------------------------------------
pdf(file.path(DIR_FIGURES, "sft_reference_unsigned.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft_ref_unsigned$fitIndices[, 1],
     -sign(sft_ref_unsigned$fitIndices[, 3]) * sft_ref_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = paste0("Scale independence (", LABEL_REFERENCE, " Unsigned)"))
text(sft_ref_unsigned$fitIndices[, 1],
     -sign(sft_ref_unsigned$fitIndices[, 3]) * sft_ref_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

plot(sft_ref_unsigned$fitIndices[, 1],
     sft_ref_unsigned$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste0("Mean connectivity (", LABEL_REFERENCE, " Unsigned)"))
text(sft_ref_unsigned$fitIndices[, 1],
     sft_ref_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

dev.off()

# ---- Test plots --------------------------------------------------------------
pdf(file.path(DIR_FIGURES, "sft_test.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft_test_signed$fitIndices[, 1],
     -sign(sft_test_signed$fitIndices[, 3]) * sft_test_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = paste0("Scale independence (", LABEL_TEST, " Signed)"))
text(sft_test_signed$fitIndices[, 1],
     -sign(sft_test_signed$fitIndices[, 3]) * sft_test_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

plot(sft_test_signed$fitIndices[, 1],
     sft_test_signed$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste0("Mean connectivity (", LABEL_TEST, " Signed)"))
text(sft_test_signed$fitIndices[, 1],
     sft_test_signed$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

dev.off()

pdf(file.path(DIR_FIGURES, "sft_test_unsigned.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft_test_unsigned$fitIndices[, 1],
     -sign(sft_test_unsigned$fitIndices[, 3]) * sft_test_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n",
     main = paste0("Scale independence (", LABEL_TEST, " Unsigned)"))
text(sft_test_unsigned$fitIndices[, 1],
     -sign(sft_test_unsigned$fitIndices[, 3]) * sft_test_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

plot(sft_test_unsigned$fitIndices[, 1],
     sft_test_unsigned$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste0("Mean connectivity (", LABEL_TEST, " Unsigned)"))
text(sft_test_unsigned$fitIndices[, 1],
     sft_test_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

dev.off()

cat("SFT diagnostic plots saved to", DIR_FIGURES, "\n")
cat("STOP HERE. Inspect the plots and fill in the POWER_* variables in Section 0.\n")


# ==============================================================================
# SECTION 10: NETWORK CONSTRUCTION
# ==============================================================================
# blockwiseModules() is the master function that executes the full pipeline:
# correlation -> adjacency (soft-thresholding) -> TOM -> hierarchical clustering
# -> dynamic tree cut -> module merging.
#
# The TOM (Topological Overlap Matrix) is saved to disk automatically so that
# if the session crashes, you do not need to recompute it.
#
# This section will run ONLY if all four power variables have been set.

if (any(is.na(c(POWER_REF_SIGNED, POWER_REF_UNSIGNED,
                POWER_TEST_SIGNED, POWER_TEST_UNSIGNED)))) {
  stop("One or more POWER_* variables are still NA. ",
       "Complete Section 8-9 and fill them in before running Section 10.")
}

cat("Building SIGNED Reference network...\n")
netwk_ref_signed <- blockwiseModules(
  t(filtered_reference),
  power           = POWER_REF_SIGNED,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = DEEP_SPLIT,
  minModuleSize   = MIN_MODULE_SIZE,
  mergeCutHeight  = MERGE_CUT_HEIGHT,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = paste0(LABEL_REFERENCE, "_signed_TOM"),
  verbose         = 3
)

cat("Building UNSIGNED Reference network...\n")
netwk_ref_unsigned <- blockwiseModules(
  t(filtered_reference),
  power           = POWER_REF_UNSIGNED,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = DEEP_SPLIT,
  minModuleSize   = MIN_MODULE_SIZE,
  mergeCutHeight  = MERGE_CUT_HEIGHT,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = paste0(LABEL_REFERENCE, "_unsigned_TOM"),
  verbose         = 3
)

cat("Building SIGNED Test network...\n")
netwk_test_signed <- blockwiseModules(
  t(filtered_test),
  power           = POWER_TEST_SIGNED,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = DEEP_SPLIT,
  minModuleSize   = MIN_MODULE_SIZE,
  mergeCutHeight  = MERGE_CUT_HEIGHT,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = paste0(LABEL_TEST, "_signed_TOM"),
  verbose         = 3
)

cat("Building UNSIGNED Test network...\n")
netwk_test_unsigned <- blockwiseModules(
  t(filtered_test),
  power           = POWER_TEST_UNSIGNED,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = DEEP_SPLIT,
  minModuleSize   = MIN_MODULE_SIZE,
  mergeCutHeight  = MERGE_CUT_HEIGHT,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = paste0(LABEL_TEST, "_unsigned_TOM"),
  verbose         = 3
)

# ---- Convert numeric module IDs to WGCNA colour strings ---------------------
colors_ref_signed    <- labels2colors(netwk_ref_signed$colors)
colors_ref_unsigned  <- labels2colors(netwk_ref_unsigned$colors)
colors_test_signed   <- labels2colors(netwk_test_signed$colors)
colors_test_unsigned <- labels2colors(netwk_test_unsigned$colors)

# ---- Module size summary tables ---------------------------------------------
cat("\n--- Module sizes ---\n")
cat(LABEL_REFERENCE, "Signed:\n");   print(table(netwk_ref_signed$colors))
cat(LABEL_REFERENCE, "Unsigned:\n"); print(table(netwk_ref_unsigned$colors))
cat(LABEL_TEST, "Signed:\n");        print(table(netwk_test_signed$colors))
cat(LABEL_TEST, "Unsigned:\n");      print(table(netwk_test_unsigned$colors))


# ==============================================================================
# SECTION 11: CLUSTER DENDROGRAMS
# ==============================================================================
# Cluster dendrograms visualise the hierarchical clustering of genes by
# topological overlap dissimilarity (1 - TOM). The colour bar underneath
# shows each gene's module assignment. Lower branch merger heights indicate
# tighter co-expression. The grey band represents unassigned (noise) genes.

cat("Generating cluster dendrograms...\n")

pdf(file.path(DIR_FIGURES, "dendrograms.pdf"), width = 14, height = 8)
par(mfrow = c(2, 2))

plotDendroAndColors(
  netwk_ref_signed$dendrograms[[1]],
  colors_ref_signed[netwk_ref_signed$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = paste0(LABEL_REFERENCE, " Network - SIGNED")
)

plotDendroAndColors(
  netwk_ref_unsigned$dendrograms[[1]],
  colors_ref_unsigned[netwk_ref_unsigned$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = paste0(LABEL_REFERENCE, " Network - UNSIGNED")
)

plotDendroAndColors(
  netwk_test_signed$dendrograms[[1]],
  colors_test_signed[netwk_test_signed$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = paste0(LABEL_TEST, " Network - SIGNED")
)

plotDendroAndColors(
  netwk_test_unsigned$dendrograms[[1]],
  colors_test_unsigned[netwk_test_unsigned$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = paste0(LABEL_TEST, " Network - UNSIGNED")
)

par(mfrow = c(1, 1))
dev.off()

cat("Dendrograms saved.\n")


# ==============================================================================
# SECTION 12: MODULE EIGENGENES AND CORRELATION HEATMAPS
# ==============================================================================
# A Module Eigengene (ME) is the first principal component of a module's gene
# expression matrix -- a single summary vector that captures the module's
# dominant expression trend across all samples. Correlating eigengenes between
# modules reveals which pathways co-activate or antagonise each other.

cat("Computing module eigengenes...\n")

compute_and_clean_MEs <- function(expression_matrix, color_vector) {
  MEs <- moduleEigengenes(t(expression_matrix), colors = color_vector)$eigengenes
  MEs <- orderMEs(MEs)
  colnames(MEs) <- gsub("ME", "", colnames(MEs))
  return(MEs)
}

MEs_ref_signed    <- compute_and_clean_MEs(filtered_reference, colors_ref_signed)
MEs_ref_unsigned  <- compute_and_clean_MEs(filtered_reference, colors_ref_unsigned)
MEs_test_signed   <- compute_and_clean_MEs(filtered_test,      colors_test_signed)
MEs_test_unsigned <- compute_and_clean_MEs(filtered_test,      colors_test_unsigned)

cor_ref_signed    <- cor(MEs_ref_signed)
cor_ref_unsigned  <- cor(MEs_ref_unsigned)
cor_test_signed   <- cor(MEs_test_signed)
cor_test_unsigned <- cor(MEs_test_unsigned)

# ---- Heatmaps ----------------------------------------------------------------
cat("Saving eigengene correlation heatmaps...\n")

library(gplots)
heatmap_palette <- colorRampPalette(c("blue", "white", "red"))(50)

pdf(file.path(DIR_FIGURES, "eigengene_heatmaps.pdf"), width = 8, height = 8)

heatmap.2(cor_ref_signed,
          main          = paste0(LABEL_REFERENCE, " Module Correlation (SIGNED)"),
          trace         = "none",
          col           = heatmap_palette,
          key           = TRUE,
          key.title     = "Correlation",
          density.info  = "none",
          denscol       = NA)

heatmap.2(cor_ref_unsigned,
          main          = paste0(LABEL_REFERENCE, " Module Correlation (UNSIGNED)"),
          trace         = "none",
          col           = heatmap_palette,
          key           = TRUE,
          key.title     = "Correlation",
          density.info  = "none",
          denscol       = NA)

heatmap.2(cor_test_signed,
          main          = paste0(LABEL_TEST, " Module Correlation (SIGNED)"),
          trace         = "none",
          col           = heatmap_palette,
          key           = TRUE,
          key.title     = "Correlation",
          density.info  = "none",
          denscol       = NA)

heatmap.2(cor_test_unsigned,
          main          = paste0(LABEL_TEST, " Module Correlation (UNSIGNED)"),
          trace         = "none",
          col           = heatmap_palette,
          key           = TRUE,
          key.title     = "Correlation",
          density.info  = "none",
          denscol       = NA)

dev.off()
cat("Heatmaps saved.\n")


# ==============================================================================
# SECTION 13: EXPORT NETWORKS TO CYTOSCAPE FORMAT
# ==============================================================================
# For each network configuration, this section writes two files:
#   Edges_<prefix>.txt : gene pairs with their TOM connection weights
#   Nodes_<prefix>.txt : gene IDs, module colours, and HGNC gene symbols
#
# The TOM is re-computed here from the expression matrix so that the edge
# weights are accurate regardless of saved TOM file format.

# EDIT: If your genes use a different ID type (e.g. SYMBOL, ENTREZID),
# change the keytype argument in mapIds() accordingly.

cat("Exporting networks to Cytoscape format...\n")

export_to_cytoscape <- function(expression_matrix,
                                power_val,
                                net_type,
                                colors_vector,
                                prefix,
                                output_dir) {
  
  cat("  Processing TOM for", prefix, "...\n")
  
  TOM_mat <- TOMsimilarityFromExpr(
    t(expression_matrix),
    power       = power_val,
    networkType = net_type
  )
  rownames(TOM_mat) <- rownames(expression_matrix)
  colnames(TOM_mat) <- rownames(expression_matrix)
  
  # Build long-format edge list and apply weight threshold
  edge_list <- as.data.frame(TOM_mat) %>%
    mutate(gene1 = rownames(.)) %>%
    pivot_longer(-gene1, names_to = "gene2", values_to = "weight") %>%
    filter(gene1 != gene2, weight > TOM_THRESHOLD) %>%
    # Remove duplicate A-B / B-A pairs
    mutate(pair = pmap_chr(list(gene1, gene2),
                           ~ paste(sort(c(...)), collapse = "_"))) %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  cat("  Translating Ensembl IDs to gene symbols for", prefix, "...\n")
  
  # EDIT: Change keytype if your IDs are not Ensembl
  edge_list$gene1.name <- mapIds(org_db,
                                 keys    = edge_list$gene1,
                                 column  = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
  edge_list$gene2.name <- mapIds(org_db,
                                 keys    = edge_list$gene2,
                                 column  = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
  
  # Fall back to Ensembl ID if no symbol found
  edge_list$gene1.name <- ifelse(is.na(edge_list$gene1.name),
                                 edge_list$gene1, edge_list$gene1.name)
  edge_list$gene2.name <- ifelse(is.na(edge_list$gene2.name),
                                 edge_list$gene2, edge_list$gene2.name)
  
  # Build node list
  node_list <- data.frame(
    GeneID      = rownames(expression_matrix),
    ModuleColor = colors_vector
  )
  node_list$GeneSymbol <- mapIds(org_db,
                                 keys    = node_list$GeneID,
                                 column  = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
  node_list$GeneSymbol <- ifelse(is.na(node_list$GeneSymbol),
                                 node_list$GeneID, node_list$GeneSymbol)
  
  # Write files
  write.table(edge_list,
              file      = file.path(output_dir, paste0("Edges_", prefix, ".txt")),
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE)
  
  write.table(node_list,
              file      = file.path(output_dir, paste0("Nodes_", prefix, ".txt")),
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE)
  
  cat("  Exported:", prefix, "\n")
}

export_to_cytoscape(filtered_reference, POWER_REF_SIGNED,
                    "signed",   colors_ref_signed,
                    paste0(LABEL_REFERENCE, "_Signed"),    DIR_MODULES)

export_to_cytoscape(filtered_reference, POWER_REF_UNSIGNED,
                    "unsigned", colors_ref_unsigned,
                    paste0(LABEL_REFERENCE, "_Unsigned"),  DIR_MODULES)

export_to_cytoscape(filtered_test, POWER_TEST_SIGNED,
                    "signed",   colors_test_signed,
                    paste0(LABEL_TEST, "_Signed"),         DIR_MODULES)

export_to_cytoscape(filtered_test, POWER_TEST_UNSIGNED,
                    "unsigned", colors_test_unsigned,
                    paste0(LABEL_TEST, "_Unsigned"),       DIR_MODULES)

cat("All Cytoscape files saved to", DIR_MODULES, "\n")


# ==============================================================================
# SECTION 14: MODULE PRESERVATION ANALYSIS
# ==============================================================================
# modulePreservation() tests whether the co-expression architecture of each
# reference (e.g. Tumour) module is maintained in the test (e.g. Normal)
# dataset. It scrambles gene labels 100 times to build a null distribution,
# then reports a Z-summary score for each module:
#   Z < 2  : no preservation (module structure collapses in test tissue)
#   2-10   : moderate preservation (core wiring survives with some changes)
#   Z > 10 : strong preservation (module is intact; likely a housekeeping pathway)
#
# The "gold" module (a random-gene control injected by WGCNA) is removed before
# plotting because its preservation score is artificial.

cat("Running module preservation analysis (this can take several hours)...\n")

# ---- Signed preservation -----------------------------------------------------
multiData_signed <- list(
  Reference = list(data = t(filtered_reference)),
  Test      = list(data = t(filtered_test))
)

ref_colors_signed  <- labels2colors(netwk_ref_signed$colors)
names(ref_colors_signed) <- names(netwk_ref_signed$colors)

test_colors_signed <- labels2colors(netwk_test_signed$colors)
names(test_colors_signed) <- names(netwk_test_signed$colors)

multiColor_signed <- list(
  Reference = ref_colors_signed,
  Test      = test_colors_signed
)

# Sanity check
stopifnot(
  all(names(multiColor_signed$Reference) %in% rownames(filtered_reference)),
  all(names(multiColor_signed$Test)      %in% rownames(filtered_test))
)
cat("Sanity checks passed for signed preservation.\n")

preservation_signed <- modulePreservation(
  multiData         = multiData_signed,
  multiColor        = multiColor_signed,
  referenceNetworks = 1,
  nPermutations     = N_PERMUTATIONS,
  randomSeed        = RANDOM_SEED,
  verbose           = 3
)

stats_signed <- preservation_signed$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test
stats_signed <- stats_signed[rownames(stats_signed) != "gold", ]

# ---- Unsigned preservation ---------------------------------------------------
multiData_unsigned <- list(
  Reference = list(data = t(filtered_reference)),
  Test      = list(data = t(filtered_test))
)

ref_colors_unsigned  <- labels2colors(netwk_ref_unsigned$colors)
names(ref_colors_unsigned) <- names(netwk_ref_unsigned$colors)

test_colors_unsigned <- labels2colors(netwk_test_unsigned$colors)
names(test_colors_unsigned) <- names(netwk_test_unsigned$colors)

multiColor_unsigned <- list(
  Reference = ref_colors_unsigned,
  Test      = test_colors_unsigned
)

stopifnot(
  all(names(multiColor_unsigned$Reference) %in% rownames(filtered_reference)),
  all(names(multiColor_unsigned$Test)      %in% rownames(filtered_test))
)
cat("Sanity checks passed for unsigned preservation.\n")

preservation_unsigned <- modulePreservation(
  multiData         = multiData_unsigned,
  multiColor        = multiColor_unsigned,
  referenceNetworks = 1,
  nPermutations     = N_PERMUTATIONS,
  randomSeed        = RANDOM_SEED,
  verbose           = 3
)

stats_unsigned <- preservation_unsigned$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test
stats_unsigned <- stats_unsigned[rownames(stats_unsigned) != "gold", ]

cat("Module preservation analysis complete.\n")

# ---- Z-summary dotplot helper ------------------------------------------------
plot_preservation <- function(stats_df, title_str, output_path) {
  plot_data <- data.frame(
    Module    = rownames(stats_df),
    Z_summary = as.numeric(stats_df$Zsummary.pres)
  )
  plot_data <- subset(plot_data, Module != "grey")
  
  p <- ggplot(plot_data, aes(x = Z_summary,
                             y = reorder(Module, Z_summary))) +
    geom_point(aes(color = Module), size = 5) +
    scale_color_identity() +
    geom_vline(xintercept = 2,  linetype = "dashed",
               color = "blue", linewidth = 0.8) +
    geom_vline(xintercept = 10, linetype = "dashed",
               color = "red",  linewidth = 0.8) +
    labs(
      title    = title_str,
      subtitle = paste0("Reference: ", LABEL_REFERENCE,
                        "  |  Test: ", LABEL_TEST),
      x = "Z-summary Score",
      y = "Module"
    ) +
    annotate("text", x = 0.8,  y = 0.6, label = "No Preservation",
             color = "blue",  angle = 90, size = 3.2) +
    annotate("text", x = 5.5,  y = 0.6, label = "Moderate",
             color = "purple", angle = 90, size = 3.2) +
    annotate("text", x = 12.5, y = 0.6, label = "Strong",
             color = "red",   angle = 90, size = 3.2) +
    theme_minimal(base_size = 12)
  
  ggsave(output_path, p, width = 6, height = 5, dpi = 300)
  return(p)
}

print(plot_preservation(
  stats_signed,
  "Module Preservation (SIGNED)",
  file.path(DIR_FIGURES, "preservation_signed.pdf")
))

print(plot_preservation(
  stats_unsigned,
  "Module Preservation (UNSIGNED)",
  file.path(DIR_FIGURES, "preservation_unsigned.pdf")
))

cat("Preservation plots saved.\n")


# ==============================================================================
# SECTION 15: EXTRACT MODULE GENE LISTS
# ==============================================================================

cat("Extracting gene lists per module...\n")

module_names_ref <- unique(ref_colors_signed)

genes_in_modules <- lapply(module_names_ref, function(mod) {
  names(ref_colors_signed[ref_colors_signed == mod])
})
names(genes_in_modules) <- module_names_ref

cat("Module gene counts (reference, signed):\n")
print(sapply(genes_in_modules, length))


# ==============================================================================
# SECTION 16: GO ENRICHMENT ANALYSIS FOR ALL MODULES
# ==============================================================================
# For each non-grey module, run enrichment against all three GO ontologies
# (Biological Process, Cellular Component, Molecular Function) and save the
# results as CSV files. Then generate an individual dotplot for each module.
#
# How to read a GO enrichment dotplot:
#   Y-axis : GO term descriptions (ordered by Gene Ratio, top = most enriched)
#   X-axis : Gene Ratio = fraction of module genes annotated to that term
#   Dot size : absolute number of module genes in the term
#   Dot colour : adjusted p-value (red = most significant, blue = less significant)
#   A term with a large, dark-red dot at a high Gene Ratio is the strongest
#   functional signature of the module.

# ---- GO enrichment helper function -------------------------------------------
run_go_enrichment <- function(gene_list, ontology, output_csv) {
  
  result <- tryCatch({
    enrichGO(
      gene          = gene_list,
      OrgDb         = org_db,
      # EDIT: Change keyType if your gene IDs are not Ensembl
      keyType       = "ENSEMBL",
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = GO_PVAL_CUTOFF,
      qvalueCutoff  = GO_QVAL_CUTOFF
    )
  }, error = function(e) {
    message("  No significant terms found for this combination.")
    return(NULL)
  })
  
  if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
    write.csv(as.data.frame(result), file = output_csv, row.names = TRUE)
  }
  
  return(result)
}

# ---- Dotplot helper ----------------------------------------------------------
save_go_dotplot <- function(go_result_csv, module_name, ontology, output_png) {
  
  if (!file.exists(go_result_csv)) return(invisible(NULL))
  go_data <- read.csv(go_result_csv)
  if (nrow(go_data) == 0) return(invisible(NULL))
  
  top_n  <- head(go_data[order(go_data$p.adjust), ], GO_TOP_CATEGORY)
  
  # Convert GeneRatio string (e.g. "12/250") to numeric
  top_n$GeneRatio_num <- sapply(top_n$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  
  p <- ggplot(top_n, aes(x = GeneRatio_num,
                         y = reorder(Description, GeneRatio_num))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title    = paste0("Top ", GO_TOP_CATEGORY, " GO-", ontology,
                        " Terms (", toupper(module_name), " Module)"),
      subtitle = paste0("Reference: ", LABEL_REFERENCE, " Network"),
      x        = "Gene Ratio",
      y        = "GO Term",
      color    = "Adjusted p-value",
      size     = "Gene Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9, face = "bold"),
      plot.title  = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(output_png, p, width = 9, height = 6, dpi = 300)
  return(p)
}

# ---- Run enrichment for all non-grey modules ---------------------------------
modules_to_test <- names(genes_in_modules)[names(genes_in_modules) != "grey"]
ontologies      <- c("BP", "CC", "MF")

cat("Running GO enrichment for", length(modules_to_test),
    "modules across 3 ontologies...\n")

for (mod in modules_to_test) {
  cat(" Module:", mod, "\n")
  
  for (ont in ontologies) {
    
    csv_path <- file.path(DIR_ENRICH,
                          paste0("GO_", ont, "_", mod, ".csv"))
    png_path <- file.path(DIR_ENRICH,
                          paste0("dotplot_", ont, "_", mod, ".png"))
    
    go_result <- run_go_enrichment(genes_in_modules[[mod]], ont, csv_path)
    
    if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
      save_go_dotplot(csv_path, mod, ont, png_path)
      cat("   GO-", ont, ": enriched terms found and plotted.\n", sep = "")
    } else {
      cat("   GO-", ont, ": no significant terms.\n", sep = "")
    }
  }
}

cat("GO enrichment complete. Results saved to", DIR_ENRICH, "\n")


# ==============================================================================
# SECTION 17: VENN DIAGRAM OF SHARED GO TERMS (two-module comparison)
# ==============================================================================
# Compare the biological process terms enriched in two modules to reveal
# whether their functional signatures are independent or overlapping.
#
# EDIT: Change MODULE_A and MODULE_B to the two modules you want to compare.
# Typically compare the two largest functional modules (e.g. turquoise and blue).

library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR)

# EDIT: Replace with the names of the two modules you want to compare.
MODULE_A <- "turquoise"
MODULE_B <- "blue"

csv_A <- file.path(DIR_ENRICH, paste0("GO_BP_", MODULE_A, ".csv"))
csv_B <- file.path(DIR_ENRICH, paste0("GO_BP_", MODULE_B, ".csv"))

if (file.exists(csv_A) && file.exists(csv_B)) {
  
  terms_A <- read.csv(csv_A)$Description
  terms_B <- read.csv(csv_B)$Description
  
  png(filename = file.path(DIR_FIGURES,
                           paste0("venn_", MODULE_A, "_vs_", MODULE_B, ".png")),
      width = 1800, height = 1800, res = 300)
  
  grid.newpage()
  venn_plot <- draw.pairwise.venn(
    area1      = length(terms_A),
    area2      = length(terms_B),
    cross.area = length(intersect(terms_A, terms_B)),
    category   = c(paste0(MODULE_A, " Module"), paste0(MODULE_B, " Module")),
    fill       = c(MODULE_A, MODULE_B),
    alpha      = c(0.45, 0.45),
    cat.pos    = c(-20, 20),
    cat.dist   = c(0.05, 0.05),
    scaled     = TRUE
  )
  grid.draw(venn_plot)
  dev.off()
  
  cat("Venn diagram saved for", MODULE_A, "vs", MODULE_B, "\n")
  cat("Shared GO-BP terms:", length(intersect(terms_A, terms_B)), "\n")
  
} else {
  cat("Venn diagram skipped: one or both module GO-BP CSV files not found.\n")
}


# ==============================================================================
# SECTION 18: HUB GENE IDENTIFICATION
# ==============================================================================
# Hub genes are the master regulators of each module: the genes with the
# highest module membership (MM), defined as the Pearson correlation between
# a gene's expression profile and the module's eigengene. Hub genes are the
# primary targets for follow-up validation experiments and literature searches.

cat("Identifying hub genes for each module...\n")

# Compute module membership (correlation of each gene with each module eigengene)
MM_ref <- cor(t(filtered_reference), MEs_ref_signed)

# Get the top N hub genes per module
N_HUBS <- 10  # EDIT: how many hub genes to report per module

hub_genes_list <- lapply(modules_to_test, function(mod) {
  
  col_name <- paste0("ME", mod)
  if (!col_name %in% colnames(MM_ref)) {
    col_name <- mod  # if colnames were cleaned
  }
  
  if (!col_name %in% colnames(MM_ref)) return(NULL)
  
  mm_scores  <- MM_ref[, col_name]
  top_genes  <- sort(abs(mm_scores), decreasing = TRUE)[1:N_HUBS]
  gene_ids   <- names(top_genes)
  
  gene_symbols <- tryCatch(
    mapIds(org_db, keys = gene_ids, column = "SYMBOL",
           keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) setNames(gene_ids, gene_ids)
  )
  
  data.frame(
    EnsemblID  = gene_ids,
    GeneSymbol = gene_symbols,
    MM_score   = top_genes,
    Module     = mod
  )
})

hub_genes_df <- do.call(rbind, hub_genes_list[!sapply(hub_genes_list, is.null)])

hub_csv_path <- file.path(DIR_ENRICH, "hub_genes_top10_per_module.csv")
write.csv(hub_genes_df, hub_csv_path, row.names = FALSE)

cat("Hub genes saved to", hub_csv_path, "\n")
print(hub_genes_df)


# ==============================================================================
# SECTION 19: SAVE WORKSPACE
# ==============================================================================
# Save the entire R session so all objects can be reloaded without re-running
# the full pipeline (which includes computationally expensive steps).

save.image(file = "WGCNA_workspace.RData")
cat("Workspace saved to WGCNA_workspace.RData\n")
cat("\nPipeline complete.\n")
cat("Key output locations:\n")
cat("  Figures   :", DIR_FIGURES, "\n")
cat("  Modules   :", DIR_MODULES, "\n")
cat("  GO results:", DIR_ENRICH,  "\n")
cat("  Workspace : WGCNA_workspace.RData\n")


# ==============================================================================
# END OF TEMPLATE
# ==============================================================================