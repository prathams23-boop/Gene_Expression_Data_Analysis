# Gene Expression Data Analysis

A complete, end-to-end pipeline implemented in R for retrieving gene expression data from The Cancer Genome Atlas (TCGA), preprocessing it, and performing Weighted Gene Co-expression Network Analysis (WGCNA). This workflow is designed to identify highly correlated gene modules and pinpoint key hub genes associated with specific clinical traits.

## 🚀 Features

- **Data Retrieval & Cleaning:** Automates the downloading, filtering, and normalization of TCGA RNA-Seq data.
- **WGCNA Implementation:** Standardized template to calculate soft-thresholding power, construct co-expression networks, and detect gene modules.
- **Module-Trait Relationships:** Correlates clinical metadata with co-expression modules.
- **Hub Gene Identification:** Extracts top network hubs based on module membership and gene significance.

------------------------------------------------------------------------

## 📁 Repository Structure

- `TCGA_Retrieval_Cleaning.R` – Script handling data extraction, sample filtering, outlier removal, and expression data normalization.
- `WGCNA_Template.R` – The core pipeline script containing network construction, module detection, and visualization code.
- `LICENSE` – MIT License.

------------------------------------------------------------------------

## 🛠️ Prerequisites & Installation

To run these scripts locally, ensure you have R (version 4.0 or higher recommended) installed along with the required Bioconductor and CRAN packages.

Open your R console and install the core dependencies:

\`\`\`R \# Install CRAN packages install.packages(c("tidyverse", "flashClust", "dynamictreeCut", "WGCNA"))

# Install Bioconductor packages (if required for data retrieval)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install(c("TCGAbiolinks", "DESeq2", "limma"))
