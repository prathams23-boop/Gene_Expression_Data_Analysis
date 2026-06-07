# TCGA-BRCA Retrieval and Cleaning

library(TCGAbiolinks)
library(SummarizedExperiment) 

TARGET_PROJECT <- "TCGA-BRCA" 

ge_query <- GDCquery(
  project       = TARGET_PROJECT,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification"
)

setwd("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01")

GDCdownload(
  ge_query,
  method         = "api",
  directory      = "C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01",
  files.per.chunk = 10
)


data_dir <- "C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01"
setwd(data_dir)


tryCatch({
  test_prepare <- GDCprepare(ge_query, directory = data_dir, summarizedExperiment = TRUE)
  message("--> Integrity Check Passed: All downloaded files are 100% healthy.")
  rm(test_prepare)
}, error = function(e) {
  message("\n[CRITICAL] Corrupted file detected: ", e$message)
  message("Delete the offending subfolder in: ", data_dir)
  stop("Pipeline halted due to data corruption.", call. = FALSE)
})


# fetching metadata
metadata_table <- getResults(ge_query)


message("  MANUAL CHECK: CRITICAL TISSUE LABELS IN THIS COHORT  ")
print(table(metadata_table$sample_type))


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

# SAFETY NET: Explicit warning if the project lacks a normal control arm 
if (length(normal_barcodes) == 0) {
  stop("\n[CRITICAL ERROR] Stop! No normal control samples were discovered using standard labels.\n",
       "Look at the printed table above. If normal tissues are labeled differently,\n",
       "add that exact spelling string to 'VALID_NORMAL_TYPES' on Line 37.")
} else {
  message("--> Success! Found ", length(normal_barcodes), " control samples and ", 
          length(tumor_barcodes), " tumor samples.")
}


# compile and extract raw expression matrix
list.files("C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01/TCGA-BRCA", 
           recursive = TRUE) |> length()

ge_data <- GDCprepare(ge_query, 
                      summarizedExperiment = TRUE,
                      directory = "C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run01" )

# Pull out raw counts (Rows = Genes, Columns = Long Tracking Barcodes)
raw_matrix <- assay(ge_data)
expression_df <- as.data.frame(raw_matrix)

# data cleaning and id harmonization
message("\n>>> STEP 6: Formatting and sanitizing sample headers...")

# Fix Suffix Mismatch: Trim everything following the initial vial code
# Example: Converts 'TCGA-CN-4723-01A-01D-...' down to 'TCGA-CN-4723-01A'
colnames(expression_df) <- substr(colnames(expression_df), 1, 16)
expression_df <- expression_df[, !duplicated(colnames(expression_df))]

# Standardize column characters (fixes R converting dashes into periods behind the scenes)
colnames(expression_df) <- gsub("\\.", "-", colnames(expression_df))
normal_barcodes         <- gsub("\\.", "-", normal_barcodes)
tumor_barcodes          <- gsub("\\.", "-", tumor_barcodes)




# Filter columns using the %in% matching criteria
normal_expression <- expression_df[, colnames(expression_df) %in% normal_barcodes, drop = FALSE]
tumor_expression  <- expression_df[, colnames(expression_df) %in% tumor_barcodes, drop = FALSE]

# final quality control checks

# Most packages (DESeq2/WGCNA) break instantly if there is an empty cell (NA)
if (any(is.na(normal_expression)) || any(is.na(tumor_expression))) {
  warning("[QC WARNING] Blank cells (NAs) detected! Consider running na.omit() on matrices.")
} else {
  message("--> QC Pass: No missing values (NAs) found in files.")
}

# Print dimensions summary

message("                 FINAL PROCESSING SUMMARY                ")

message("  Normal Data Matrix: ", nrow(normal_expression), " Genes across ", ncol(normal_expression), " Samples.")
message("  Tumor Data Matrix:  ", nrow(tumor_expression), " Genes across ", ncol(tumor_expression), " Samples.")
message("=========================================================")

# export clean spreadsheets


write.csv(normal_expression, paste0(TARGET_PROJECT, "_clean_expression_normal.csv"), row.names = TRUE)
write.csv(tumor_expression,  paste0(TARGET_PROJECT, "_clean_expression_tumor.csv"),  row.names = TRUE)

message("\n[SUCCESS] Pipeline complete! Your data frames are formatted and ready for DESeq2 processing and WGCNA.")
# ==============================================================================