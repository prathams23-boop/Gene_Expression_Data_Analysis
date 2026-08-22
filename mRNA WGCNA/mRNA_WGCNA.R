
# WGCNA pipeline for TCGA-BRCA mRNA expression data
# tumour-only cohort (963 samples), PAM50 subtype used as sample trait
# gene selection is a MAD based variance filter


library(readxl)
library(DESeq2)
library(WGCNA)
library(matrixStats)
library(gplots)
library(ggplot2)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()


MRNA_XLSX      <- "/Users/PrathamShah/Desktop/IIITH/Codes/WGCNA_Pipeline/Gene_Expression_Data_Analysis/mRNA_filter1.xlsx"
OUT_DIR        <- "/Users/PrathamShah/Desktop/IIITH/Codes/WGCNA_Pipeline/Gene_Expression_Data_Analysis/mRNA WGCNA"

# low-expression filter: keep a gene if it has at least MIN_COUNT reads
# in at least MIN_FRAC of samples
MIN_COUNT      <- 10
MIN_FRAC       <- 0.50

# number of genes carried into WGCNA after ranking by MAD
# set to NA to keep everything that survives the expression filter
N_GENES_WGCNA  <- NA

# single-block network construction, do not leave this at the default 5000
MAX_BLOCK_SIZE <- 12000

N_HUBS         <- 20

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),    showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "enrich"),    showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "modules"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "hub_genes"), showWarnings = FALSE)


# loading data


mrna_raw <- as.data.frame(read_excel(MRNA_XLSX, sheet = 1))

rownames(mrna_raw) <- mrna_raw$Sample_ID
mrna_raw$Sample_ID <- NULL

cat("Raw dimensions (samples x columns):", dim(mrna_raw), "\n")


stopifnot("Subtype" %in% colnames(mrna_raw))

subtype_mrna <- mrna_raw$Subtype
names(subtype_mrna) <- rownames(mrna_raw)

# NOTE: "TNBC" in this file is the Basal PAM50 class under a different
# name. Relabeling so downstream tables match paper_BRCA_Subtype_PAM50.
subtype_mrna[subtype_mrna == "TNBC"] <- "Basal"
print(table(subtype_mrna))

expr_mrna <- mrna_raw[, setdiff(colnames(mrna_raw), "Subtype"), drop = FALSE]
expr_mrna <- as.matrix(expr_mrna)

# typecasting for DESeq2
mode(expr_mrna) <- "numeric"
stopifnot(all(expr_mrna == round(expr_mrna)))
expr_mrna <- round(expr_mrna)
storage.mode(expr_mrna) <- "integer"

cat("Count matrix (samples x genes):", dim(expr_mrna), "\n")
cat("Fraction of zero entries:", round(mean(expr_mrna == 0), 4), "\n")
cat("Library size median / min / max:",
    median(rowSums(expr_mrna)), min(rowSums(expr_mrna)), max(rowSums(expr_mrna)), "\n")


counts_mrna <- t(expr_mrna)   # genes x samples
cat("DESeq2 input (genes x samples):", dim(counts_mrna), "\n")


# low expression filter

keep_expressed <- rowMeans(counts_mrna >= MIN_COUNT) >= MIN_FRAC
cat("Genes passing expression filter:", sum(keep_expressed),
    "of", nrow(counts_mrna), "\n")

counts_mrna <- counts_mrna[keep_expressed, ]


# normalisation and vst
# design = ~1 because there is no group comparison here. DESeq2 is used
# only for size factor estimation and the dispersion fit that VST needs.

coldata_mrna <- data.frame(
  row.names = colnames(counts_mrna),
  Subtype   = factor(subtype_mrna[colnames(counts_mrna)])
)
stopifnot(identical(rownames(coldata_mrna), colnames(counts_mrna)))

dds_mrna <- DESeqDataSetFromMatrix(
  countData = counts_mrna,
  colData   = coldata_mrna,
  design    = ~ 1
)

dds_mrna <- estimateSizeFactors(dds_mrna)

cat("Size factor range:", round(range(sizeFactors(dds_mrna)), 3), "\n")
summary(sizeFactors(dds_mrna))

# saving the size factors, the lncRNA script will reuse them, since the
# lncRNA submatrix is too sparse to estimate depth reliably on its own
saveRDS(sizeFactors(dds_mrna),
        file.path(OUT_DIR, "mrna_size_factors.rds"))

dds_mrna <- estimateDispersions(dds_mrna)


vsd_mrna  <- varianceStabilizingTransformation(dds_mrna, blind = TRUE)
vst_mrna  <- assay(vsd_mrna)   # genes x samples

cat("VST matrix (genes x samples):", dim(vst_mrna), "\n")


saveRDS(vst_mrna, file.path(OUT_DIR, "mrna_vst_all_expressed.rds"))


# GENE SELECTION BY MAD
gene_mad <- rowMads(vst_mrna)
names(gene_mad) <- rownames(vst_mrna)

png(file.path(OUT_DIR, "figures", "mrna_mad_distribution.png"),
    width = 1400, height = 900, res = 150)
hist(gene_mad, breaks = 80,
     main = "mRNA gene-level MAD (VST scale)",
     xlab = "Median absolute deviation")
if (!is.na(N_GENES_WGCNA)) {
  abline(v = sort(gene_mad, decreasing = TRUE)[min(N_GENES_WGCNA, length(gene_mad))],
         col = "red", lty = 2, lwd = 2)
}
dev.off()

if (is.na(N_GENES_WGCNA)) {
  selected_genes <- names(gene_mad)
} else {
  selected_genes <- names(sort(gene_mad, decreasing = TRUE))[1:min(N_GENES_WGCNA, length(gene_mad))]
}

filtered_mrna <- vst_mrna[selected_genes, ]
cat("Genes going into WGCNA:", nrow(filtered_mrna), "\n")

# WGCNA wants samples x genes. Keep this orientation from here on.
datExpr_mrna <- t(filtered_mrna)
cat("WGCNA input (samples x genes):", dim(datExpr_mrna), "\n")





gsg_mrna <- goodSamplesGenes(datExpr_mrna, verbose = 3)
if (!gsg_mrna$allOK) {
  cat("Removing", sum(!gsg_mrna$goodGenes),   "genes and",
      sum(!gsg_mrna$goodSamples), "samples flagged by goodSamplesGenes\n")
  datExpr_mrna <- datExpr_mrna[gsg_mrna$goodSamples, gsg_mrna$goodGenes]
}

# sample clustering to spot gross outlier arrays
sampleTree_mrna <- hclust(dist(datExpr_mrna), method = "average")

png(file.path(OUT_DIR, "figures", "mrna_sample_clustering.png"),
    width = 2400, height = 900, res = 150)
plot(sampleTree_mrna,
     main = "mRNA sample clustering (outlier detection)",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4,
     labels = FALSE)
# draw a candidate cut line, adjust after looking at the plot
CUT_HEIGHT_MRNA <- NA  
if (!is.na(CUT_HEIGHT_MRNA)) abline(h = CUT_HEIGHT_MRNA, col = "red", lty = 2)
dev.off()

if (!is.na(CUT_HEIGHT_MRNA)) {
  clust_mrna   <- cutreeStatic(sampleTree_mrna, cutHeight = CUT_HEIGHT_MRNA, minSize = 10)
  keep_samples <- clust_mrna == names(sort(table(clust_mrna[clust_mrna != 0]), decreasing = TRUE))[1]
  cat("Dropping", sum(!keep_samples), "outlier samples\n")
  datExpr_mrna <- datExpr_mrna[keep_samples, ]
}

subtype_mrna <- subtype_mrna[rownames(datExpr_mrna)]
cat("Final cohort:", nrow(datExpr_mrna), "samples,", ncol(datExpr_mrna), "genes\n")


# sft selection

powers <- c(1:10, seq(from = 12, to = 20, by = 2))

sft_mrna_signed <- pickSoftThreshold(
  datExpr_mrna,
  powerVector = powers,
  networkType = "signed",
  verbose     = 5
)

sft_mrna_unsigned <- pickSoftThreshold(
  datExpr_mrna,
  powerVector = powers,
  networkType = "unsigned",
  verbose     = 5
)

print(sft_mrna_signed$fitIndices)
print(sft_mrna_unsigned$fitIndices)

# scale independence and mean connectivity plots
cex1 <- 0.9

png(file.path(OUT_DIR, "figures", "mrna_soft_threshold.png"),
    width = 2000, height = 1600, res = 160)
par(mfrow = c(2, 2))

plot(sft_mrna_unsigned$fitIndices[, 1],
     -sign(sft_mrna_unsigned$fitIndices[, 3]) * sft_mrna_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence (mRNA, Unsigned)")
text(sft_mrna_unsigned$fitIndices[, 1],
     -sign(sft_mrna_unsigned$fitIndices[, 3]) * sft_mrna_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red", lty = 2)
abline(h = 0.80, col = "blue", lty = 3)

plot(sft_mrna_signed$fitIndices[, 1],
     -sign(sft_mrna_signed$fitIndices[, 3]) * sft_mrna_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence (mRNA, Signed)")
text(sft_mrna_signed$fitIndices[, 1],
     -sign(sft_mrna_signed$fitIndices[, 3]) * sft_mrna_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red", lty = 2)
abline(h = 0.80, col = "blue", lty = 3)

plot(sft_mrna_unsigned$fitIndices[, 1], sft_mrna_unsigned$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity (mRNA, Unsigned)")
text(sft_mrna_unsigned$fitIndices[, 1], sft_mrna_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

plot(sft_mrna_signed$fitIndices[, 1], sft_mrna_signed$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity (mRNA, Signed)")
text(sft_mrna_signed$fitIndices[, 1], sft_mrna_signed$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

par(mfrow = c(1, 1))
dev.off()

# helper function: lowest power crossing an R^2 target with negative slope
pick_power <- function(sft, target = 0.90, fallback) {
  fi  <- sft$fitIndices
  r2  <- -sign(fi[, 3]) * fi[, 2]
  ok  <- which(r2 >= target)
  if (length(ok) == 0) {
    cat("  no power reaches R^2 =", target, ", falling back to", fallback, "\n")
    return(fallback)
  }
  fi[min(ok), 1]
}

power_mrna_signed   <- pick_power(sft_mrna_signed,   0.90, fallback = 12)
power_mrna_unsigned <- pick_power(sft_mrna_unsigned, 0.90, fallback = 6)

cat("Chosen power, signed:",   power_mrna_signed,   "\n")
cat("Chosen power, unsigned:", power_mrna_unsigned, "\n")

# for manual selection
# power_mrna_signed   <- 9
# power_mrna_unsigned <- 4

# network construction

cor <- WGCNA::cor 

# first run
netwk_mrna_signed <- blockwiseModules(
  datExpr_mrna,
  power           = power_mrna_signed,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = 2,
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "mrna_signed_TOM"),
  verbose         = 3
)

netwk_mrna_unsigned <- blockwiseModules(
  datExpr_mrna,
  power           = power_mrna_unsigned,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = 2,
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "mrna_unsigned_TOM"),
  verbose         = 3
)

cat("Signed, first run module sizes:\n");   print(table(netwk_mrna_signed$colors))
cat("Unsigned, first run module sizes:\n"); print(table(netwk_mrna_unsigned$colors))

# tuned run
# balanced modules. Target is roughly 5 to 12 non-grey modules.
netwk_mrna_signed_1 <- blockwiseModules(
  datExpr_mrna,
  power           = power_mrna_signed,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = 1,
  minModuleSize   = 50,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "mrna_signed_1_TOM"),
  verbose         = 3
)

netwk_mrna_unsigned_1 <- blockwiseModules(
  datExpr_mrna,
  power           = power_mrna_unsigned,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = 1,
  minModuleSize   = 50,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "mrna_unsigned_1_TOM"),
  verbose         = 3
)

# confirm single-block construction
cat("Signed blocks:",   length(netwk_mrna_signed_1$blockGenes),
    "| Unsigned blocks:", length(netwk_mrna_unsigned_1$blockGenes), "\n")


colors_mrna_signed   <- labels2colors(netwk_mrna_signed_1$colors)
colors_mrna_unsigned <- labels2colors(netwk_mrna_unsigned_1$colors)

names(colors_mrna_signed)   <- colnames(datExpr_mrna)
names(colors_mrna_unsigned) <- colnames(datExpr_mrna)

cat("Signed, tuned module sizes:\n");   print(table(colors_mrna_signed))
cat("Unsigned, tuned module sizes:\n"); print(table(colors_mrna_unsigned))

write.csv(as.data.frame(table(colors_mrna_signed)),
          file.path(OUT_DIR, "tables", "mrna_module_sizes_signed.csv"), row.names = FALSE)
write.csv(as.data.frame(table(colors_mrna_unsigned)),
          file.path(OUT_DIR, "tables", "mrna_module_sizes_unsigned.csv"), row.names = FALSE)

# gene to module assignment table
write.csv(
  data.frame(
    GeneID        = colnames(datExpr_mrna),
    ModuleSigned  = colors_mrna_signed,
    ModuleUnsigned = colors_mrna_unsigned
  ),
  file.path(OUT_DIR, "tables", "mrna_gene_module_assignment.csv"), row.names = FALSE
)


# cluster dendrograms

png(file.path(OUT_DIR, "figures", "mrna_dendrograms.png"),
    width = 2400, height = 1600, res = 160)
par(mfrow = c(2, 1))

plotDendroAndColors(
  netwk_mrna_signed_1$dendrograms[[1]],
  colors_mrna_signed[netwk_mrna_signed_1$blockGenes[[1]]],
  "Module colours",
  main = paste0("mRNA Network - SIGNED (tuned, power ", power_mrna_signed, ")"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

plotDendroAndColors(
  netwk_mrna_unsigned_1$dendrograms[[1]],
  colors_mrna_unsigned[netwk_mrna_unsigned_1$blockGenes[[1]]],
  "Module colours",
  main = paste0("mRNA Network - UNSIGNED (tuned, power ", power_mrna_unsigned, ")"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

par(mfrow = c(1, 1))
dev.off()


# module eigengenes and correlation heatmaps

MEs_mrna_signed <- moduleEigengenes(datExpr_mrna, colors = colors_mrna_signed)$eigengenes
MEs_mrna_signed <- MEs_mrna_signed[, colnames(MEs_mrna_signed) != "MEgrey", drop = FALSE]
MEs_mrna_signed <- orderMEs(MEs_mrna_signed)
colnames(MEs_mrna_signed) <- gsub("^ME", "", colnames(MEs_mrna_signed))

MEs_mrna_unsigned <- moduleEigengenes(datExpr_mrna, colors = colors_mrna_unsigned)$eigengenes
MEs_mrna_unsigned <- MEs_mrna_unsigned[, colnames(MEs_mrna_unsigned) != "MEgrey", drop = FALSE]
MEs_mrna_unsigned <- orderMEs(MEs_mrna_unsigned)
colnames(MEs_mrna_unsigned) <- gsub("^ME", "", colnames(MEs_mrna_unsigned))

cor_mrna_signed   <- cor(MEs_mrna_signed)
cor_mrna_unsigned <- cor(MEs_mrna_unsigned)

png(file.path(OUT_DIR, "figures", "mrna_eigengene_heatmaps.png"),
    width = 2000, height = 1000, res = 150)
par(mfrow = c(1, 2))

heatmap.2(cor_mrna_signed,
          main = "mRNA Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

heatmap.2(cor_mrna_unsigned,
          main = "mRNA Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

par(mfrow = c(1, 1))
dev.off()

saveRDS(MEs_mrna_signed,   file.path(OUT_DIR, "mrna_MEs_signed.rds"))
saveRDS(MEs_mrna_unsigned, file.path(OUT_DIR, "mrna_MEs_unsigned.rds"))


# module trait relationship (trial, open for comments)

subtype_factor <- factor(subtype_mrna[rownames(datExpr_mrna)])
trait_mat <- model.matrix(~ 0 + subtype_factor)
colnames(trait_mat) <- levels(subtype_factor)
rownames(trait_mat) <- rownames(datExpr_mrna)


cat("Trait matrix:\n"); print(colSums(trait_mat))

module_trait_analysis <- function(MEs, traits, tag) {
  
  nSamples <- nrow(MEs)
  mtCor    <- cor(MEs, traits, use = "pairwise.complete.obs")
  mtP      <- corPvalueStudent(mtCor, nSamples)
  
  textMatrix <- paste0(signif(mtCor, 2), "\n(", signif(mtP, 1), ")")
  dim(textMatrix) <- dim(mtCor)
  
  png(file.path(OUT_DIR, "figures", paste0("mrna_module_trait_", tag, ".png")),
      width = 1600, height = 1400, res = 150)
  par(mar = c(6, 9, 3, 3))
  labeledHeatmap(
    Matrix        = mtCor,
    xLabels       = colnames(traits),
    yLabels       = paste0("ME", rownames(mtCor)),
    ySymbols      = rownames(mtCor),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = textMatrix,
    setStdMargins = FALSE,
    cex.text      = 0.7,
    zlim          = c(-1, 1),
    main          = paste0("Module-Subtype Relationships (mRNA, ", tag, ")")
  )
  dev.off()
  
  out <- as.data.frame(mtCor)
  colnames(out) <- paste0("cor_", colnames(out))
  outP <- as.data.frame(mtP)
  colnames(outP) <- paste0("p_", colnames(outP))
  res <- cbind(Module = rownames(mtCor), out, outP)
  
  write.csv(res, file.path(OUT_DIR, "tables",
                           paste0("mrna_module_trait_", tag, ".csv")),
            row.names = FALSE)
  res
}

mt_signed   <- module_trait_analysis(MEs_mrna_signed,   trait_mat, "signed")
mt_unsigned <- module_trait_analysis(MEs_mrna_unsigned, trait_mat, "unsigned")

print(mt_signed)


# go enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

dir.create(file.path(OUT_DIR, "enrich", "GO_Signed"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "enrich", "GO_Unsigned"), showWarnings = FALSE)

# Ensembl version suffixes must be stripped, org.Hs.eg.db only accepts
# bare IDs and fails silently otherwise
perform_go_enrichment <- function(gene_list, universe_list, ontology, output_path) {
  
  gene_clean     <- gsub("\\..*", "", gene_list)
  universe_clean <- gsub("\\..*", "", universe_list)
  
  go_results <- enrichGO(
    gene          = gene_clean,
    universe      = universe_clean,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    write.csv(as.data.frame(go_results), file = output_path, row.names = TRUE)
  }
  go_results
}

plot_go_for_module <- function(go_bp, go_cc, go_mf, mod_name, network_type, out_dir) {
  
  plots <- list()
  for (nm in c("BP", "CC", "MF")) {
    obj <- switch(nm, BP = go_bp, CC = go_cc, MF = go_mf)
    lab <- switch(nm,
                  BP = "GO Biological Process",
                  CC = "GO Cellular Component",
                  MF = "GO Molecular Function")
    if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
      plots[[nm]] <- dotplot(obj, showCategory = 10, font.size = 9, label_format = 60) +
        scale_size_continuous(range = c(2, 7)) +
        theme_minimal() +
        ggtitle(paste0(lab, " (", mod_name, " - ", network_type, ")"))
    }
  }
  
  if (length(plots) > 0) {
    combined <- ggarrange(plotlist = plots, ncol = 1, nrow = length(plots))
    out_file <- file.path(out_dir, paste0("GO_combined_", mod_name, ".png"))
    ggsave(out_file, combined, width = 10, height = 5 * length(plots), dpi = 300)
    cat("  saved dotplot for", mod_name, "\n")
  } else {
    cat("  no significant GO terms in any ontology for module:", mod_name, "\n")
  }
}

run_go_for_network <- function(colors_vec, network_type, out_dir) {
  
  gene_universe <- names(colors_vec)
  mods <- setdiff(unique(colors_vec), "grey")
  cat("\nModules to test (", network_type, "):", paste(mods, collapse = ", "), "\n")
  
  results <- list()
  for (mod in mods) {
    cat("\nProcessing module:", mod, "(", network_type, ")\n")
    gene_list <- names(colors_vec)[colors_vec == mod]
    results[[mod]] <- list()
    
    for (ont in c("BP", "CC", "MF")) {
      out_path <- file.path(out_dir, paste0("GO_", ont, "_", mod, ".csv"))
      results[[mod]][[ont]] <- tryCatch(
        perform_go_enrichment(gene_list, gene_universe, ont, out_path),
        error = function(e) { cat("  ", ont, "failed:", e$message, "\n"); NULL }
      )
    }
    
    plot_go_for_module(results[[mod]][["BP"]], results[[mod]][["CC"]],
                       results[[mod]][["MF"]], mod, network_type, out_dir)
  }
  results
}

go_mrna_signed   <- run_go_for_network(colors_mrna_signed,   "Signed",
                                       file.path(OUT_DIR, "enrich", "GO_Signed"))
go_mrna_unsigned <- run_go_for_network(colors_mrna_unsigned, "Unsigned",
                                       file.path(OUT_DIR, "enrich", "GO_Unsigned"))


# hub gene identification

extract_hub_genes <- function(datExpr, MEs, colors_vec, tag) {
  
  MM <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
  cat("\n", tag, "module membership matrix:", dim(MM), "\n")
  
  hub_list <- list()
  
  for (mod in colnames(MM)) {
    
    kme_scores <- MM[[mod]]
    names(kme_scores) <- colnames(datExpr)
    
    # restrict to assigned members
    mod_genes  <- names(colors_vec)[colors_vec == mod]
    kme_scores <- kme_scores[mod_genes]
    
    n_take <- min(N_HUBS, length(kme_scores))
    if (n_take < N_HUBS) {
      cat("  NOTE: module", mod, "has only", length(kme_scores), "genes\n")
    }
    
    # rank by signed kME, true members are positive
    top_ids  <- names(sort(kme_scores, decreasing = TRUE))[1:n_take]
    top_kmes <- kme_scores[top_ids]
    
    top_clean <- gsub("\\..*", "", top_ids)
    gene_symbols <- tryCatch(
      mapIds(org.Hs.eg.db, keys = top_clean, column = "SYMBOL",
             keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) setNames(top_ids, top_ids)
    )
    gene_symbols[is.na(gene_symbols)] <- top_ids[is.na(gene_symbols)]
    
    hub_list[[mod]] <- data.frame(
      Module     = mod,
      EnsemblID  = top_ids,
      GeneSymbol = as.character(gene_symbols),
      kME        = round(as.numeric(top_kmes), 4),
      row.names  = NULL
    )
    
    cat("  Module", mod, "top hub:", as.character(gene_symbols)[1],
        "(kME =", round(top_kmes[1], 3), ")\n")
  }
  
  hubs <- do.call(rbind, hub_list)
  rownames(hubs) <- NULL
  
  write.csv(hubs, file.path(OUT_DIR, "hub_genes",
                            paste0("mrna_hub_genes_", tag, "_top", N_HUBS, ".csv")),
            row.names = FALSE)
  
  # no gene should now appear in two modules
  dup_check <- hubs %>%
    group_by(GeneSymbol) %>%
    summarise(n_modules = n_distinct(Module), .groups = "drop") %>%
    filter(n_modules > 1)
  cat("Genes in more than one", tag, "module:", nrow(dup_check), "(should be 0)\n")
  if (nrow(dup_check) > 0) print(dup_check)
  
  hubs
}

hub_mrna_signed   <- extract_hub_genes(datExpr_mrna, MEs_mrna_signed,
                                       colors_mrna_signed,   "signed")
hub_mrna_unsigned <- extract_hub_genes(datExpr_mrna, MEs_mrna_unsigned,
                                       colors_mrna_unsigned, "unsigned")

# kME bar charts per module
plot_kme_bars <- function(hubs, tag) {
  p <- ggplot(hubs, aes(x = reorder(GeneSymbol, kME), y = kME, fill = Module)) +
    geom_col() +
    scale_fill_manual(values = setNames(unique(hubs$Module), unique(hubs$Module))) +
    coord_flip() +
    facet_wrap(~ Module, scales = "free_y") +
    labs(title = paste0("Top ", N_HUBS, " hub genes by kME (mRNA, ", tag, ")"),
         x = NULL, y = "kME") +
    theme_bw() + theme(legend.position = "none", axis.text.y = element_text(size = 6))
  
  ggsave(file.path(OUT_DIR, "figures", paste0("mrna_kme_bars_", tag, ".png")),
         p, width = 14, height = 10, dpi = 300)
}

plot_kme_bars(hub_mrna_signed,   "signed")
plot_kme_bars(hub_mrna_unsigned, "unsigned")



# save session objects

cor <- stats::cor   

save(datExpr_mrna, subtype_mrna,
     power_mrna_signed, power_mrna_unsigned,
     netwk_mrna_signed_1, netwk_mrna_unsigned_1,
     colors_mrna_signed, colors_mrna_unsigned,
     MEs_mrna_signed, MEs_mrna_unsigned,
     mt_signed, mt_unsigned,
     hub_mrna_signed, hub_mrna_unsigned,
     file = file.path(OUT_DIR, "mrna_wgcna_session.RData"))

cat("\nmRNA WGCNA pipeline complete. Outputs in:", OUT_DIR, "\n")
sessionInfo()