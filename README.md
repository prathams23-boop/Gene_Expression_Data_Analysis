# Gene Expression Data Analysis

A complete, end-to-end pipeline implemented in R for retrieving gene expression data from The Cancer Genome Atlas (TCGA), preprocessing it, and performing Weighted Gene Co-expression Network Analysis (WGCNA). This workflow is designed to identify highly correlated gene modules and pinpoint key hub genes associated with specific clinical traits.

This repository is a product of my work carried out as a summer intern (since May 2026) at the International Institute of Information Technology (IIIT), Hyderabad under the supervision of Professor Nita Parekh.

## 🚀 Features

- **Data Retrieval & Cleaning:** Automates the downloading, filtering, and normalization of TCGA RNA-Seq data.
- **WGCNA Implementation:** Standardized template to calculate soft-thresholding power, construct co-expression networks, and detect gene modules.
- **Module-Trait Relationships:** Correlates clinical metadata with co-expression modules.
- **Hub Gene Identification:** Extracts top network hubs based on module membership and gene significance.

------------------------------------------------------------------------

## 📁 Repository Structure

- `TCGA_BRCA_Retrieval.R` – Script handling data extraction, sample filtering, outlier removal, and expression data normalization.
- `TCGA_BRCA_Full_Pipeline.R` – The core pipeline script containing network construction, module detection, and visualization code.
- `WGCNA_Tutorial_Downloads.R` – All the important packages you need for this
- `WGCNA_Protocol.pdf` – The protocol that I referred for the same
- `LICENSE` – MIT License.

------------------------------------------------------------------------

## 🛠️ Prerequisites & Notes

To run these scripts locally, ensure you have R (version 4.0 or higher recommended) installed along with the required Bioconductor and CRAN packages.

Please refer to WGCNA_Tutorial_Downloads.R for instructions on installing all the required packages.

The protocol I have attached (which I referred to while learning DESeq2 and WGCNA) does not include analyses for both signed and unsigned networks. My version, however, includes a detailed analysis of both these networks. While picking the soft threshold, I have considered multiple candidate powers for both the networks and have constructed multiple networks accordingly. The optimum power and network was then chosen on the basis of module separation as visible in the dendrogram plots.

Another change from the protocol is that they have used OSCC data while I have worked here on BRCA data.

In the detailed report that I have attached, you will find comprehensive explanations of the theory behind the whole code and algorithm when and where needed, along with proper explanations on how to interpret key figures.

## 🛠️ How to use this?

Using these scripts is easy, in the retrieval script, just change the project name according to the TCGA project data you will be working with. In the full pipeline script, you just need to change some variables and commands initially which take the number of features of the data as input, rest of it will be pretty much the same. Going through the entire thing here will surely be helpful before moving on to using it for your project.

#### **Thank you, and good luck!**
