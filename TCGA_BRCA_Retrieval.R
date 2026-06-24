# tcga-brca full retrieval and cleaning pipeline

library(TCGAbiolinks)
library(SummarizedExperiment)

TARGET_PROJECT <- "TCGA-BRCA"
DATA_DIR       <- "C:/Users/Pratham Shah/Downloads/IIITH/TCGA_BRCA_Data_Run02"

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(DATA_DIR)


# query and download

ge_query <- GDCquery(
  project               = TARGET_PROJECT,
  data.category         = "Transcriptome Profiling",
  data.type             = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type         = "STAR - Counts",
  platform              = "Illumina",
  data.format           = "tsv",
  access                = "open"
)

GDCdownload(
  ge_query,
  method          = "api",
  directory       = DATA_DIR,
  files.per.chunk = 5
)


# just a check to see if files are corrupted, only need to use if download was interrupted 
# and had to be started again, otherwise its fine

tryCatch({
  test_prepare <- GDCprepare(ge_query, directory = DATA_DIR, summarizedExperiment = TRUE)
  message("--> Integrity check passed: all downloaded files are healthy.")
  rm(test_prepare)
}, error = function(e) {
  message("\n[CRITICAL] Corrupted file detected: ", e$message)
  stop("Pipeline halted due to data corruption.", call. = FALSE)
})


# extract gene annotation

ge_data <- GDCprepare(ge_query,
                      summarizedExperiment = TRUE,
                      directory = DATA_DIR)

# Pull gene annotation from rowRanges metadata
gene_annotation <- data.frame(
  gene_id   = ge_data@rowRanges@elementMetadata@listData$gene_id,
  gene_name = ge_data@rowRanges@elementMetadata$gene_name,
  gene_type = ge_data@rowRanges@elementMetadata$gene_type,
  stringsAsFactors = FALSE
)

write.csv(gene_annotation, "Gene_Annotation_All.csv", row.names = FALSE)
message("--> Saved: Gene_Annotation_All.csv  (", nrow(gene_annotation), " genes)")

# Keep only protein-coding genes and remove duplicate gene names
protein_coding <- gene_annotation[gene_annotation$gene_type == "protein_coding", ]
protein_coding <- protein_coding[!duplicated(protein_coding$gene_name), ]

write.csv(protein_coding, "Gene_Annotation_ProteinCoding.csv", row.names = FALSE)
message("--> Saved: Gene_Annotation_ProteinCoding.csv  (",
        nrow(protein_coding), " protein-coding genes)")



# BUILD RAW COUNTS MATRIX (protein-coding genes only)
#            First two columns: ensembl_id, gene_name
#            Remaining columns: one per sample (16-char barcode)

raw_matrix <- assay(ge_data)

# Subset rows to protein-coding Ensembl IDs
brca_expr <- as.data.frame(raw_matrix[protein_coding$gene_id, ])

# Trim sample barcodes to 16-character form (e.g. TCGA-A2-A0CM-01A)
colnames(brca_expr) <- substr(colnames(brca_expr), 1, 16)

# Drop any duplicate sample columns that arose from barcode trimming
brca_expr <- brca_expr[, !duplicated(colnames(brca_expr))]

# Standardise dashes (R sometimes converts them to dots internally)
colnames(brca_expr) <- gsub("\\.", "-", colnames(brca_expr))

# Prepend Ensembl ID and gene name as identifier columns
brca_expr <- cbind(
  ensembl_id = protein_coding$gene_id,
  gene_name  = protein_coding$gene_name,
  brca_expr
)

write.csv(brca_expr, "TCGA-BRCA_raw_counts_protein_coding_all_samples.csv",
          row.names = FALSE)
message("--> Saved: TCGA-BRCA_raw_counts_protein_coding_all_samples.csv")
message("    Dimensions: ", nrow(brca_expr), " genes x ",
        ncol(brca_expr) - 2, " samples")


# CLINICAL METADATA EXTRACTION
#            (a) Receptor status (ER / PR / HER2) via BCR Biotab
#            (b) PAM50 subtype + sample type via colDataPrepare


# BCR Biotab download (receptor status)
sample_query <- GDCquery(
  project       = TARGET_PROJECT,
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)

GDCdownload(sample_query, method = "api", directory = DATA_DIR)

clinical_brca  <- GDCprepare(sample_query, directory = DATA_DIR)
patient_biotab <- clinical_brca$clinical_patient_brca

# Rows 1-2 are header artefacts in BCR Biotab exports -- always skip them
patient_biotab <- patient_biotab[3:nrow(patient_biotab), ]

patient_receptor <- patient_biotab[, c("bcr_patient_barcode",
                                       "er_status_by_ihc",
                                       "pr_status_by_ihc",
                                       "her2_status_by_ihc")]
colnames(patient_receptor) <- c("patient", "ER_Status", "PR_Status", "HER2_Status")

# additional clinical fields (age, menopause status, tumor stage, race)
# patient_biotab already has the header-artefact rows removed above.
# Run print(colnames(patient_biotab)) once and confirm these exact names exist
# in your export before trusting this block; BCR Biotab column names can drift
# slightly across cohorts/versions (e.g. ajcc_pathologic_tumor_stage is
# sometimes named tumor_stage instead).
print(colnames(patient_biotab))

patient_extra <- patient_biotab[, c("bcr_patient_barcode",
                                    "age_at_diagnosis",
                                    "menopause_status",
                                    "ajcc_pathologic_tumor_stage",
                                    "race")]
colnames(patient_extra) <- c("patient", "Age", "Menopause_Status",
                             "Tumor_Stage", "Race")

patient_extra$Age <- as.numeric(patient_extra$Age)

# colDataPrepare for PAM50 subtype + sample_type
# Use the 16-char sample IDs already in the expression matrix
sample_ids <- colnames(brca_expr)[-(1:2)]   # drop the two annotation columns

clinical_col <- colDataPrepare(sample_ids)
clinical_col <- clinical_col[, c("patient", "sample", "sample_type",
                                 "paper_BRCA_Subtype_PAM50")]

# Merge both clinical tables on patient barcode 
clinical_merged <- merge(patient_receptor, clinical_col, by = "patient")
clinical_merged <- merge(clinical_merged, patient_extra, by = "patient")  # NEW

write.csv(clinical_merged, "TCGA-BRCA_clinical_metadata.csv", row.names = FALSE)
message("--> Saved: TCGA-BRCA_clinical_metadata.csv  (",
        nrow(clinical_merged), " samples)")



# REMOVE NORMAL-LIKE SAMPLES, THEN SPLIT TUMOR vs NORMAL


message("\n  PAM50 subtype distribution before filtering:")
print(table(clinical_merged$paper_BRCA_Subtype_PAM50, useNA = "ifany"))
message("\n  Sample type distribution before filtering:")
print(table(clinical_merged$sample_type, useNA = "ifany"))

# Identify sample-level barcodes of PAM50 Normal-Like tumour samples
normallike_samples <- clinical_merged$sample[
  !is.na(clinical_merged$paper_BRCA_Subtype_PAM50) &
    clinical_merged$paper_BRCA_Subtype_PAM50 == "Normal"
]

message("\n--> Removing ", length(normallike_samples),
        " Normal-Like PAM50 samples from expression matrix.")

# Build the filtered expression matrix
annotation_cols  <- c("ensembl_id", "gene_name")
expr_sample_cols <- colnames(brca_expr)[-(1:2)]
keep_samples     <- expr_sample_cols[!expr_sample_cols %in% normallike_samples]

brca_expr_filtered <- brca_expr[, c(annotation_cols, keep_samples)]

# Split using sample_type from clinical metadata
VALID_TUMOR_TYPES  <- c("Primary Tumor")
VALID_NORMAL_TYPES <- c("Solid Tissue Normal",
                        "Blood Derived Normal",
                        "Peripheral Blood Lymphocyte")

tumor_ids  <- clinical_merged$sample[
  clinical_merged$sample_type %in% VALID_TUMOR_TYPES &
    !clinical_merged$sample %in% normallike_samples
]

normal_ids <- clinical_merged$sample[
  clinical_merged$sample_type %in% VALID_NORMAL_TYPES
]

tumor_expr  <- brca_expr_filtered[, c(annotation_cols,
                                      intersect(keep_samples, tumor_ids))]
normal_expr <- brca_expr_filtered[, c(annotation_cols,
                                      intersect(keep_samples, normal_ids))]


# final quality checks and export


for (mat_name in c("tumor_expr", "normal_expr")) {
  mat       <- get(mat_name)
  expr_only <- mat[, -(1:2)]
  if (any(is.na(expr_only))) {
    warning("[QC WARNING] NAs detected in ", mat_name)
  } else {
    message("--> QC pass: no NAs in ", mat_name)
  }
}

message("\n  FINAL SUMMARY")
message("  Protein-coding genes:            ", nrow(tumor_expr))
message("  Tumor samples (Normal-Like removed): ", ncol(tumor_expr) - 2)
message("  Normal tissue samples:           ", ncol(normal_expr) - 2)


write.csv(tumor_expr,  "TCGA-BRCA_protein_coding_tumor.csv",  row.names = FALSE)
write.csv(normal_expr, "TCGA-BRCA_protein_coding_normal.csv", row.names = FALSE)

message("\n[SUCCESS] All outputs written to: ", DATA_DIR)
message("  Gene_Annotation_All.csv")
message("  Gene_Annotation_ProteinCoding.csv")
message("  TCGA-BRCA_raw_counts_protein_coding_all_samples.csv")
message("  TCGA-BRCA_clinical_metadata.csv")
message("  TCGA-BRCA_protein_coding_tumor.csv")
message("  TCGA-BRCA_protein_coding_normal.csv")
