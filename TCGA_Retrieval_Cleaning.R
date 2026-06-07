# ==============================================================================
# MASTER TEMPLATE: TCGA RNA-Seq Data Retrieval & Universal Cleaning Pipeline
# Description: Downloads raw counts from GDC, aligns with metadata, cleans IDs,
#              handles variant tissue labels, and exports isolated CSV matrices.
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Load Required Libraries
# ------------------------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment) # Needed to extract data layers smoothly

# ------------------------------------------------------------------------------
# STEP 2: Configure Project Settings
# ------------------------------------------------------------------------------
# Change this string to study different cancers! Examples:
# "TCGA-HNSC" (Oral/Head & Neck), "TCGA-BRCA" (Breast), "TCGA-LUAD" (Lung)
TARGET_PROJECT <- "TCGA-HNSC" 

# ------------------------------------------------------------------------------
# STEP 3: Query and Download from GDC
# ------------------------------------------------------------------------------
message("\n>>> STEP 3: Querying GDC database for ", TARGET_PROJECT, "...")

ge_query <- GDCquery(
  project       = TARGET_PROJECT,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification"
)

# Using 'api' method because it bypasses broken internal gdc-client download URLs
message(">>> Downloading files from GDC via API... (Please be patient)")
GDCdownload(ge_query, method = "api")

# ------------------------------------------------------------------------------
# STEP 4: Fetch Metadata & Handle Tissue Labels (CRITICAL ATTENTION ZONE)
# ------------------------------------------------------------------------------
message("\n>>> STEP 4: Extracting sample indexing metadata...")
metadata_table <- getResults(ge_query)

# --- THE RADAR SCREEN: Look at your console output when this runs! ---
message("=========================================================")
message("  MANUAL CHECK: CRITICAL TISSUE LABELS IN THIS COHORT  ")
message("=========================================================")
print(table(metadata_table$sample_type))
message("=========================================================")

# Define default TCGA groups
VALID_TUMOR_TYPES  <- c("Primary Tumor")
# Added alternative labels for blood/rare tissue cohorts (e.g., Leukemia)
VALID_NORMAL_TYPES <- c("Solid Tissue Normal", "Blood Derived Normal", "Peripheral Blood Lymphocyte")

# Isolate the sample barcode tracking strings
normal_barcodes <- metadata_table$sample.submitter_id[
  metadata_table$sample_type %in% VALID_NORMAL_TYPES
]

tumor_barcodes <- metadata_table$sample.submitter_id[
  metadata_table$sample_type %in% VALID_TUMOR_TYPES
]

# --- SAFETY NET: Explicit warning if the project lacks a normal control arm ---
if (length(normal_barcodes) == 0) {
  stop("\n[CRITICAL ERROR] Stop! No normal control samples were discovered using standard labels.\n",
       "Look at the printed table above. If normal tissues are labeled differently,\n",
       "add that exact spelling string to 'VALID_NORMAL_TYPES' on Line 37.")
} else {
  message("--> Success! Found ", length(normal_barcodes), " control samples and ", 
          length(tumor_barcodes), " tumor samples.")
}

# ------------------------------------------------------------------------------
# STEP 5: Compile and Extract Raw Expression Matrix
# ------------------------------------------------------------------------------
message("\n>>> STEP 5: Compiling downloaded files into an R DataFrame...")
ge_data <- GDCprepare(ge_query, summarizedExperiment = TRUE)

# Pull out raw counts (Rows = Genes, Columns = Long Tracking Barcodes)
raw_matrix <- assay(ge_data)
expression_df <- as.data.frame(raw_matrix)

# ------------------------------------------------------------------------------
# STEP 6: Data Cleaning & ID Harmonization
# ------------------------------------------------------------------------------
message("\n>>> STEP 6: Formatting and sanitizing sample headers...")

# Fix Suffix Mismatch: Trim everything following the initial vial code
# Example: Converts 'TCGA-CN-4723-01A-01D-...' down to 'TCGA-CN-4723-01A'
colnames(expression_df) <- gsub("-.*$", "", colnames(expression_df))

# Standardize column characters (fixes R converting dashes into periods behind the scenes)
colnames(expression_df) <- gsub("\\.", "-", colnames(expression_df))
normal_barcodes         <- gsub("\\.", "-", normal_barcodes)
tumor_barcodes          <- gsub("\\.", "-", tumor_barcodes)

# ------------------------------------------------------------------------------
# STEP 7: Subset Master Dataset into Isolated Dataframes
# ------------------------------------------------------------------------------
message("\n>>> STEP 7: Slicing data matrix into separate tissue environments...")

# Filter columns using the %in% matching criteria
normal_expression <- expression_df[, colnames(expression_df) %in% normal_barcodes, drop = FALSE]
tumor_expression  <- expression_df[, colnames(expression_df) %in% tumor_barcodes, drop = FALSE]

# ------------------------------------------------------------------------------
# STEP 8: Final Quality Control Check (NA Screening)
# ------------------------------------------------------------------------------
message("\n>>> STEP 8: Running final execution quality checks...")

# Most packages (DESeq2/WGCNA) break instantly if there is an empty cell (NA)
if (any(is.na(normal_expression)) || any(is.na(tumor_expression))) {
  warning("[QC WARNING] Blank cells (NAs) detected! Consider running na.omit() on matrices.")
} else {
  message("--> QC Pass: No missing values (NAs) found in files.")
}

# Print dimensions summary
message("\n=========================================================")
message("                 FINAL PROCESSING SUMMARY                ")
message("=========================================================")
message("  Normal Data Matrix: ", nrow(normal_expression), " Genes across ", ncol(normal_expression), " Samples.")
message("  Tumor Data Matrix:  ", nrow(tumor_expression), " Genes across ", ncol(tumor_expression), " Samples.")
message("=========================================================")

# ------------------------------------------------------------------------------
# STEP 9: Export Clean Spreadsheets
# ------------------------------------------------------------------------------
message("\n>>> STEP 9: Outputting clean CSV data matrices to your workspace...")

write.csv(normal_expression, paste0(TARGET_PROJECT, "_clean_expression_normal.csv"), row.names = TRUE)
write.csv(tumor_expression,  paste0(TARGET_PROJECT, "_clean_expression_tumor.csv"),  row.names = TRUE)

message("\n[SUCCESS] Pipeline complete! Your data frames are formatted and ready for Step 3 DESeq2 processing.")
# ==============================================================================