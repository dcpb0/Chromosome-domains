# Set personal library
dir.create("~/R/library", recursive = TRUE, showWarnings = FALSE)
.libPaths(c("~/R/library", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript r_deseq.r <counts_tsv> <metadata_tsv> [out_results_tsv] [out_shrunk_tsv]")
}
counts_path <- args[[1]]
metadata_path <- args[[2]]
out_results_path <- if (length(args) >= 3) args[[3]] else "DEA_results_270126.tsv"
out_shrunk_path <- if (length(args) >= 4) args[[4]] else "DEA_shrunk_results_270126.tsv"

# Install core dependencies
install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("GenomicRanges")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")

library(DESeq2)

# Load the count data and metadata
counts <- read.csv(counts_path, sep='\t', row.names=1)
metadata <- read.csv(metadata_path, sep='\t', row.names=1)
colnames(counts) <- rownames(metadata)

# Create a DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ cell_range)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Extract results and mean expression
results_list <- list()
means_list <- list()
baseline <- "0-1"
means <- as.data.frame(rowMeans(counts(dds, normalized=TRUE)[, metadata$cell_range == baseline]))
colnames(means) <- paste(baseline, "meanExpr", sep="_")
means_list[[baseline]] <- means

for (cell_range in setdiff(unique(metadata$cell_range), baseline)) {
  # Results for log2FC and adjusted p-values
  res <- results(dds, contrast=c("cell_range", cell_range, baseline))
  results_list[[cell_range]] <- as.data.frame(res)[, c("log2FoldChange", "padj")]
  colnames(results_list[[cell_range]]) <- c(paste(cell_range, "log2FoldChange", sep="_"), paste(cell_range, "padj", sep="_"))
  
  # Mean expressions for each group
  means <- as.data.frame(rowMeans(assays(dds)[["counts"]][, metadata$cell_range == cell_range]))
  colnames(means) <- paste(cell_range, "meanExpr", sep="_")
  means_list[[cell_range]] <- means
  print(colnames(results_list[[cell_range]]))
  
}

# Step 2: Concatenate all dataframes within results_list and means_list separately
results_combined <- do.call(cbind, results_list)
colnames(results_combined) <- gsub("(.*)\\..*_", "\\1_", colnames(results_combined))

means_combined <- do.call(cbind, means_list)
print(colnames(results_combined))

# Step 3: Merge the concatenated results and means dataframes
combined_res <- cbind(results_combined, means_combined)

# Save the result as a CSV file
# write.table(combined_res, file=out_results_path, sep='\t')


summary(results(dds, contrast=c("cell_range", "5-8", "0-1")))

plotDispEsts(dds)


res <- results(dds, contrast=c("cell_range", "5-8", "0-1"))
plotMA(res, ylim=c(-5,5))


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="cell_range")


res_5_8 <- lfcShrink(dds,
                     coef="cell_range_5.8_vs_0.1",
                     type="apeglm")
summary(res_5_8)

# MA plot
plotMA(res_5_8, ylim=c(-5,5), main="MA plot: 5-8 vs 0-1")


dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ cell_range)
dds <- DESeq(dds)

results_list <- list()
means_list <- list()
baseline <- "0-1"

means <- as.data.frame(rowMeans(counts(dds, normalized=TRUE)[, metadata$cell_range == baseline]))
colnames(means) <- paste(baseline, "meanExpr", sep="_")
means_list[[baseline]] <- means

for (cell_range in setdiff(unique(metadata$cell_range), baseline)) {
  
  coef_name <- paste0("cell_range_", gsub("-", ".", gsub("\\+", ".", cell_range)), "_vs_0.1")
  
  res_shrunk <- lfcShrink(dds, coef=coef_name, type="apeglm")
  
  # Extract log2FC and padj
  results_list[[cell_range]] <- as.data.frame(res_shrunk)[, c("log2FoldChange", "padj")]
  colnames(results_list[[cell_range]]) <- c(paste(cell_range, "log2FoldChange", sep="_"),
                                            paste(cell_range, "padj", sep="_"))
  
  means <- as.data.frame(rowMeans(counts(dds, normalized=TRUE)[, metadata$cell_range == cell_range]))
  colnames(means) <- paste(cell_range, "meanExpr", sep="_")
  means_list[[cell_range]] <- means
}

results_combined <- do.call(cbind, results_list)
colnames(results_combined) <- gsub("(.*)\\..*_", "\\1_", colnames(results_combined))

means_combined <- do.call(cbind, means_list)

# Merge and save
combined_res <- cbind(results_combined, means_combined)
# write.table(combined_res,
#             file=out_shrunk_path,
#             sep='\t',
#             quote=FALSE,
#             row.names=TRUE)


