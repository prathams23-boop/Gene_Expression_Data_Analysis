# WGCNA tutorial package install

library(DESeq2)
library(tidyverse)

install.packages(c(
  "WGCNA",
  "dendextend",
  "ggplots",
  "VennDiagram",
  "GO.db"
))

BiocManager::install(c(
  "DESeq2",
  "genefilter",
  "clusterProfiler",
  "org.Hs.eg.db",
  "TCGAbiolinks"
))
