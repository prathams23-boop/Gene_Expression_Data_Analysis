# WGCNA pipeline for TCGA-BRCA lncRNA expression data
# tumour-only cohort (963 samples, identical cohort to the mRNA run),
# PAM50 subtype used as sample trait
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


LNCRNA_XLSX         <- "/Users/PrathamShah/Desktop/IIITH/Codes/WGCNA_Pipeline/Gene_Expression_Data_Analysis/lncRNA_filter1.xlsx"
OUT_DIR             <- "/Users/PrathamShah/Desktop/IIITH/Codes/WGCNA_Pipeline/Gene_Expression_Data_Analysis/lncRNA WGCNA"

# size factors exported by the mRNA pipeline
MRNA_SIZE_FACTORS_RDS <- "/Users/PrathamShah/Desktop/IIITH/Codes/WGCNA_Pipeline/Gene_Expression_Data_Analysis/mRNA WGCNA/R_data_mRNA/mrna_size_factors.rds"

# low-expression filter
MIN_COUNT      <- 5
MIN_FRAC       <- 0.50

# number of genes carried into WGCNA after ranking by MAD
# setting to NA for this first run to keep everything that survives the expression filter
N_GENES_WGCNA  <- NA

# single-block network construction
MAX_BLOCK_SIZE <- 12000

N_HUBS         <- 20

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),    showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "enrich"),    showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "modules"),   showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "hub_genes"), showWarnings = FALSE)


# loading data


lncrna_raw <- as.data.frame(read_excel(LNCRNA_XLSX, sheet = 1))

rownames(lncrna_raw) <- lncrna_raw$Sample_ID
lncrna_raw$Sample_ID <- NULL

cat("Raw dimensions (samples x columns):", dim(lncrna_raw), "\n")


stopifnot("Subtype" %in% colnames(lncrna_raw))

subtype_lncrna <- lncrna_raw$Subtype
names(subtype_lncrna) <- rownames(lncrna_raw)


subtype_lncrna[subtype_lncrna == "TNBC"] <- "Basal"
print(table(subtype_lncrna))

expr_lncrna <- lncrna_raw[, setdiff(colnames(lncrna_raw), "Subtype"), drop = FALSE]
expr_lncrna <- as.matrix(expr_lncrna)

# typecasting for DESeq2
mode(expr_lncrna) <- "numeric"
stopifnot(all(expr_lncrna == round(expr_lncrna)))
expr_lncrna <- round(expr_lncrna)
storage.mode(expr_lncrna) <- "integer"

cat("Count matrix (samples x genes):", dim(expr_lncrna), "\n")
cat("Fraction of zero entries:", round(mean(expr_lncrna == 0), 4), "\n")
cat("Library size median / min / max:",
    median(rowSums(expr_lncrna)), min(rowSums(expr_lncrna)), max(rowSums(expr_lncrna)), "\n")


counts_lncrna <- t(expr_lncrna)   # genes x samples
cat("DESeq2 input (genes x samples):", dim(counts_lncrna), "\n")


# low expression filter

keep_expressed <- rowMeans(counts_lncrna >= MIN_COUNT) >= MIN_FRAC
cat("Genes passing expression filter:", sum(keep_expressed),
    "of", nrow(counts_lncrna), "\n")

counts_lncrna <- counts_lncrna[keep_expressed, ]


# normalisation and vst


coldata_lncrna <- data.frame(
  row.names = colnames(counts_lncrna),
  Subtype   = factor(subtype_lncrna[colnames(counts_lncrna)])
)
stopifnot(identical(rownames(coldata_lncrna), colnames(counts_lncrna)))

dds_lncrna <- DESeqDataSetFromMatrix(
  countData = counts_lncrna,
  colData   = coldata_lncrna,
  design    = ~ 1
)

# import mRNA size factors instead of re-estimating 
mrna_size_factors <- readRDS(MRNA_SIZE_FACTORS_RDS)

common_samples <- intersect(names(mrna_size_factors), colnames(dds_lncrna))
cat("Samples shared between mRNA size factors and lncRNA cohort:",
    length(common_samples), "of", ncol(dds_lncrna), "\n")

if (length(common_samples) < ncol(dds_lncrna)) {
  warning("Not every lncRNA sample has a matching mRNA size factor; ",
          "dropping unmatched samples so the two cohorts stay aligned.")
  dds_lncrna     <- dds_lncrna[, common_samples]
  subtype_lncrna <- subtype_lncrna[common_samples]
}

sizeFactors(dds_lncrna) <- mrna_size_factors[colnames(dds_lncrna)]

cat("Imported size factor range:", round(range(sizeFactors(dds_lncrna)), 3), "\n")
summary(sizeFactors(dds_lncrna))


dds_lncrna <- estimateDispersions(dds_lncrna)


vsd_lncrna  <- varianceStabilizingTransformation(dds_lncrna, blind = TRUE)
vst_lncrna  <- assay(vsd_lncrna)   # genes x samples

cat("VST matrix (genes x samples):", dim(vst_lncrna), "\n")


saveRDS(vst_lncrna, file.path(OUT_DIR, "lncrna_vst_all_expressed.rds"))


# GENE SELECTION BY MAD
gene_mad <- rowMads(vst_lncrna)
names(gene_mad) <- rownames(vst_lncrna)

png(file.path(OUT_DIR, "figures", "lncrna_mad_distribution.png"),
    width = 1400, height = 900, res = 150)
hist(gene_mad, breaks = 80,
     main = "lncRNA gene-level MAD (VST scale)",
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

filtered_lncrna <- vst_lncrna[selected_genes, ]
cat("Genes going into WGCNA:", nrow(filtered_lncrna), "\n")

# WGCNA wants samples x genes. Keep this orientation from here on.
datExpr_lncrna <- t(filtered_lncrna)
cat("WGCNA input (samples x genes):", dim(datExpr_lncrna), "\n")




gsg_lncrna <- goodSamplesGenes(datExpr_lncrna, verbose = 3)
if (!gsg_lncrna$allOK) {
  cat("Removing", sum(!gsg_lncrna$goodGenes),   "genes and",
      sum(!gsg_lncrna$goodSamples), "samples flagged by goodSamplesGenes\n")
  datExpr_lncrna <- datExpr_lncrna[gsg_lncrna$goodSamples, gsg_lncrna$goodGenes]
}

# sample clustering to spot gross outlier arrays
sampleTree_lncrna <- hclust(dist(datExpr_lncrna), method = "average")

png(file.path(OUT_DIR, "figures", "lncrna_sample_clustering.png"),
    width = 2400, height = 900, res = 150)
plot(sampleTree_lncrna,
     main = "lncRNA sample clustering (outlier detection)",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4,
     labels = FALSE)
# draw a candidate cut line, adjust after looking at the plot
CUT_HEIGHT_LNCRNA <- NA  # set to a numeric height to actually cut
if (!is.na(CUT_HEIGHT_LNCRNA)) abline(h = CUT_HEIGHT_LNCRNA, col = "red", lty = 2)
dev.off()

if (!is.na(CUT_HEIGHT_LNCRNA)) {
  clust_lncrna <- cutreeStatic(sampleTree_lncrna, cutHeight = CUT_HEIGHT_LNCRNA, minSize = 10)
  keep_samples <- clust_lncrna == names(sort(table(clust_lncrna[clust_lncrna != 0]), decreasing = TRUE))[1]
  cat("Dropping", sum(!keep_samples), "outlier samples\n")
  datExpr_lncrna <- datExpr_lncrna[keep_samples, ]
}

subtype_lncrna <- subtype_lncrna[rownames(datExpr_lncrna)]
cat("Final cohort:", nrow(datExpr_lncrna), "samples,", ncol(datExpr_lncrna), "genes\n")


# sft selection

powers <- c(1:10, seq(from = 12, to = 20, by = 2))

sft_lncrna_signed <- pickSoftThreshold(
  datExpr_lncrna,
  powerVector = powers,
  networkType = "signed",
  verbose     = 5
)

sft_lncrna_unsigned <- pickSoftThreshold(
  datExpr_lncrna,
  powerVector = powers,
  networkType = "unsigned",
  verbose     = 5
)

print(sft_lncrna_signed$fitIndices)
print(sft_lncrna_unsigned$fitIndices)

# scale independence and mean connectivity plots
cex1 <- 0.9

png(file.path(OUT_DIR, "figures", "lncrna_soft_threshold.png"),
    width = 2000, height = 1600, res = 160)
par(mfrow = c(2, 2))

plot(sft_lncrna_unsigned$fitIndices[, 1],
     -sign(sft_lncrna_unsigned$fitIndices[, 3]) * sft_lncrna_unsigned$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence (lncRNA, Unsigned)")
text(sft_lncrna_unsigned$fitIndices[, 1],
     -sign(sft_lncrna_unsigned$fitIndices[, 3]) * sft_lncrna_unsigned$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red", lty = 2)
abline(h = 0.80, col = "blue", lty = 3)

plot(sft_lncrna_signed$fitIndices[, 1],
     -sign(sft_lncrna_signed$fitIndices[, 3]) * sft_lncrna_signed$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence (lncRNA, Signed)")
text(sft_lncrna_signed$fitIndices[, 1],
     -sign(sft_lncrna_signed$fitIndices[, 3]) * sft_lncrna_signed$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red", lty = 2)
abline(h = 0.80, col = "blue", lty = 3)

plot(sft_lncrna_unsigned$fitIndices[, 1], sft_lncrna_unsigned$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity (lncRNA, Unsigned)")
text(sft_lncrna_unsigned$fitIndices[, 1], sft_lncrna_unsigned$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

plot(sft_lncrna_signed$fitIndices[, 1], sft_lncrna_signed$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity (lncRNA, Signed)")
text(sft_lncrna_signed$fitIndices[, 1], sft_lncrna_signed$fitIndices[, 5],
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

power_lncrna_signed   <- pick_power(sft_lncrna_signed,   0.90, fallback = 14)
power_lncrna_unsigned <- pick_power(sft_lncrna_unsigned, 0.90, fallback = 7)

cat("Chosen power, signed:",   power_lncrna_signed,   "\n")
cat("Chosen power, unsigned:", power_lncrna_unsigned, "\n")

# for manual selection
# power_lncrna_signed   <- 9
# power_lncrna_unsigned <- 4

# network construction

cor <- WGCNA::cor

# first run
netwk_lncrna_signed <- blockwiseModules(
  datExpr_lncrna,
  power           = power_lncrna_signed,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = 2,
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "lncrna_signed_TOM"),
  verbose         = 3
)

netwk_lncrna_unsigned <- blockwiseModules(
  datExpr_lncrna,
  power           = power_lncrna_unsigned,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = 2,
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "lncrna_unsigned_TOM"),
  verbose         = 3
)

cat("Signed, first run module sizes:\n");   print(table(netwk_lncrna_signed$colors))
cat("Unsigned, first run module sizes:\n"); print(table(netwk_lncrna_unsigned$colors))

# tuned run

netwk_lncrna_signed_1 <- blockwiseModules(
  datExpr_lncrna,
  power           = power_lncrna_signed,
  networkType     = "signed",
  TOMType         = "signed",
  deepSplit       = 1,
  minModuleSize   = 50,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "lncrna_signed_1_TOM"),
  verbose         = 3
)

netwk_lncrna_unsigned_1 <- blockwiseModules(
  datExpr_lncrna,
  power           = power_lncrna_unsigned,
  networkType     = "unsigned",
  TOMType         = "unsigned",
  deepSplit       = 1,
  minModuleSize   = 50,
  mergeCutHeight  = 0.25,
  maxBlockSize    = MAX_BLOCK_SIZE,
  numericLabels   = TRUE,
  saveTOMs        = TRUE,
  saveTOMFileBase = file.path(OUT_DIR, "lncrna_unsigned_1_TOM"),
  verbose         = 3
)

# confirm single-block construction
cat("Signed blocks:",   length(netwk_lncrna_signed_1$blockGenes),
    "| Unsigned blocks:", length(netwk_lncrna_unsigned_1$blockGenes), "\n")


colors_lncrna_signed   <- labels2colors(netwk_lncrna_signed$colors)
colors_lncrna_unsigned <- labels2colors(netwk_lncrna_unsigned$colors)

names(colors_lncrna_signed)   <- colnames(datExpr_lncrna)
names(colors_lncrna_unsigned) <- colnames(datExpr_lncrna)

cat("Signed, tuned module sizes:\n");   print(table(colors_lncrna_signed))
cat("Unsigned, tuned module sizes:\n"); print(table(colors_lncrna_unsigned))

write.csv(as.data.frame(table(colors_lncrna_signed)),
          file.path(OUT_DIR, "tables", "lncrna_module_sizes_signed.csv"), row.names = FALSE)
write.csv(as.data.frame(table(colors_lncrna_unsigned)),
          file.path(OUT_DIR, "tables", "lncrna_module_sizes_unsigned.csv"), row.names = FALSE)

# gene to module assignment table
write.csv(
  data.frame(
    GeneID         = colnames(datExpr_lncrna),
    ModuleSigned   = colors_lncrna_signed,
    ModuleUnsigned = colors_lncrna_unsigned
  ),
  file.path(OUT_DIR, "tables", "lncrna_gene_module_assignment.csv"), row.names = FALSE
)


# cluster dendrograms

png(file.path(OUT_DIR, "figures", "lncrna_dendrograms.png"),
    width = 2400, height = 1600, res = 160)
par(mfrow = c(2, 1))

plotDendroAndColors(
  netwk_lncrna_signed_1$dendrograms[[1]],
  colors_lncrna_signed[netwk_lncrna_signed$blockGenes[[1]]],
  "Module colours",
  main = paste0("lncRNA Network - SIGNED (tuned, power ", power_lncrna_signed, ")"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

plotDendroAndColors(
  netwk_lncrna_unsigned_1$dendrograms[[1]],
  colors_lncrna_unsigned[netwk_lncrna_unsigned$blockGenes[[1]]],
  "Module colours",
  main = paste0("lncRNA Network - UNSIGNED (tuned, power ", power_lncrna_unsigned, ")"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

par(mfrow = c(1, 1))
dev.off()


# module eigengenes and correlation heatmaps

MEs_lncrna_signed <- moduleEigengenes(datExpr_lncrna, colors = colors_lncrna_signed)$eigengenes
MEs_lncrna_signed <- MEs_lncrna_signed[, colnames(MEs_lncrna_signed) != "MEgrey", drop = FALSE]
MEs_lncrna_signed <- orderMEs(MEs_lncrna_signed)
colnames(MEs_lncrna_signed) <- gsub("^ME", "", colnames(MEs_lncrna_signed))

MEs_lncrna_unsigned <- moduleEigengenes(datExpr_lncrna, colors = colors_lncrna_unsigned)$eigengenes
MEs_lncrna_unsigned <- MEs_lncrna_unsigned[, colnames(MEs_lncrna_unsigned) != "MEgrey", drop = FALSE]
MEs_lncrna_unsigned <- orderMEs(MEs_lncrna_unsigned)
colnames(MEs_lncrna_unsigned) <- gsub("^ME", "", colnames(MEs_lncrna_unsigned))

cor_lncrna_signed   <- cor(MEs_lncrna_signed)
cor_lncrna_unsigned <- cor(MEs_lncrna_unsigned)

png(file.path(OUT_DIR, "figures", "lncrna_eigengene_heatmaps.png"),
    width = 2000, height = 1000, res = 150)
par(mfrow = c(1, 2))

heatmap.2(cor_lncrna_signed,
          main = "lncRNA Module Correlation\n(SIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

heatmap.2(cor_lncrna_unsigned,
          main = "lncRNA Module Correlation\n(UNSIGNED)",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(50),
          key = TRUE, key.title = "Correlation", key.xlab = "Value",
          density.info = "none", denscol = NA)

par(mfrow = c(1, 1))
dev.off()

saveRDS(MEs_lncrna_signed,   file.path(OUT_DIR, "lncrna_MEs_signed.rds"))
saveRDS(MEs_lncrna_unsigned, file.path(OUT_DIR, "lncrna_MEs_unsigned.rds"))


# module trait relationship (trial, open for comments)

subtype_factor <- factor(subtype_lncrna[rownames(datExpr_lncrna)])
trait_mat <- model.matrix(~ 0 + subtype_factor)
colnames(trait_mat) <- levels(subtype_factor)
rownames(trait_mat) <- rownames(datExpr_lncrna)


cat("Trait matrix:\n"); print(colSums(trait_mat))

module_trait_analysis <- function(MEs, traits, tag) {
  
  nSamples <- nrow(MEs)
  mtCor    <- cor(MEs, traits, use = "pairwise.complete.obs")
  mtP      <- corPvalueStudent(mtCor, nSamples)
  
  textMatrix <- paste0(signif(mtCor, 2), "\n(", signif(mtP, 1), ")")
  dim(textMatrix) <- dim(mtCor)
  
  png(file.path(OUT_DIR, "figures", paste0("lncrna_module_trait_", tag, ".png")),
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
    main          = paste0("Module-Subtype Relationships (lncRNA, ", tag, ")")
  )
  dev.off()
  
  out <- as.data.frame(mtCor)
  colnames(out) <- paste0("cor_", colnames(out))
  outP <- as.data.frame(mtP)
  colnames(outP) <- paste0("p_", colnames(outP))
  res <- cbind(Module = rownames(mtCor), out, outP)
  
  write.csv(res, file.path(OUT_DIR, "tables",
                           paste0("lncrna_module_trait_", tag, ".csv")),
            row.names = FALSE)
  res
}

mt_lncrna_signed   <- module_trait_analysis(MEs_lncrna_signed,   trait_mat, "signed")
mt_lncrna_unsigned <- module_trait_analysis(MEs_lncrna_unsigned, trait_mat, "unsigned")

print(mt_lncrna_signed)





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
    
    hub_list[[mod]] <- data.frame(
      Module     = mod,
      GeneID     = top_ids,
      GeneSymbol = top_ids,
      kME        = round(as.numeric(top_kmes), 4),
      row.names  = NULL
    )
    
    cat("  Module", mod, "top hub:", top_ids[1],
        "(kME =", round(top_kmes[1], 3), ")\n")
  }
  
  hubs <- do.call(rbind, hub_list)
  rownames(hubs) <- NULL
  
  write.csv(hubs, file.path(OUT_DIR, "hub_genes",
                            paste0("lncrna_hub_genes_", tag, "_top", N_HUBS, ".csv")),
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

hub_lncrna_signed   <- extract_hub_genes(datExpr_lncrna, MEs_lncrna_signed,
                                         colors_lncrna_signed,   "signed")
hub_lncrna_unsigned <- extract_hub_genes(datExpr_lncrna, MEs_lncrna_unsigned,
                                         colors_lncrna_unsigned, "unsigned")

# kME bar charts per module
plot_kme_bars <- function(hubs, tag) {
  p <- ggplot(hubs, aes(x = reorder(GeneSymbol, kME), y = kME, fill = Module)) +
    geom_col() +
    scale_fill_manual(values = setNames(unique(hubs$Module), unique(hubs$Module))) +
    coord_flip() +
    facet_wrap(~ Module, scales = "free_y") +
    labs(title = paste0("Top ", N_HUBS, " hub genes by kME (lncRNA, ", tag, ")"),
         x = NULL, y = "kME") +
    theme_bw() + theme(legend.position = "none", axis.text.y = element_text(size = 6))
  
  ggsave(file.path(OUT_DIR, "figures", paste0("lncrna_kme_bars_", tag, ".png")),
         p, width = 14, height = 10, dpi = 300)
}

plot_kme_bars(hub_lncrna_signed,   "signed")
plot_kme_bars(hub_lncrna_unsigned, "unsigned")



# save session objects

cor <- stats::cor

save(datExpr_lncrna, subtype_lncrna,
     power_lncrna_signed, power_lncrna_unsigned,
     netwk_lncrna_signed_1, netwk_lncrna_unsigned_1,
     colors_lncrna_signed, colors_lncrna_unsigned,
     MEs_lncrna_signed, MEs_lncrna_unsigned,
     mt_lncrna_signed, mt_lncrna_unsigned,
     hub_lncrna_signed, hub_lncrna_unsigned,
     file = file.path(OUT_DIR, "lncrna_wgcna_session.RData"))

cat("\nlncRNA WGCNA pipeline complete. Outputs in:", OUT_DIR, "\n")
sessionInfo()