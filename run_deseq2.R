#!/usr/bin/env Rscript
###############################################################################
# Script Name:        run_deseq2.R
# Author:             Israel Barrantes
# Email:              israel.barrantes@uni-rostock.de
#
# Version:            v2.1.0 (2026-10-07)
#
# Description:
#   Differential expression analysis (cancer vs normal) with DESeq2. Maps
#   transcript IDs to gene symbols using a GTF file, writes the DEG table and
#   generates a PCA plot and a heatmap of the top 20 genes by adjusted p-value.
#
# Required R Version:
#   R >= <version>
#
# Required Packages:
#   DESeq2       (>= <version>)
#   optparse     (>= <version>)
#
# Usage:
#   Rscript run_deseq2.R -m metadata.tsv -e counts.tsv -g refs/22.gtf \
#     [-d degs.tsv] [-p pca_plot.png,heatmap.png]
#
# Inputs:
#   --metadata / -m     Sample metadata (TSV); requires columns "sample" and
#                       "condition" (levels "normal" and "cancer")
#   --expression / -e   Count matrix (TSV); first column = transcript/gene ID,
#                       remaining columns = samples
#   --gtf / -g          GTF annotation used to map transcript IDs to gene symbols
#
# Outputs:
#   --output-degs / -d  TSV with transcript_id, gene_symbol, log2FoldChange,
#                       padj (NA padj removed, sorted by padj)
#                       [default: degs.tsv]
#   --output-plots / -p PCA plot and heatmap (PNG), comma-separated
#                       [default: pca_plot.png,heatmap.png]
#
# Notes:
#   - Genes with total counts < 10 are removed before fitting
#   - Design: ~ condition; reference level "normal"; contrast cancer vs normal
#   - Symbols not found in the GTF fall back to the unversioned transcript ID
#   - DEG counts reported on console use padj < 0.05 and |log2FC| > 1;
#     the output table itself is not filtered by these thresholds
#   - Heatmap: top 20 genes by padj, VST counts, row-wise z-scores
#
# Changelog:
#   2026-10-07 — v2.1.0 — Standard header
#   2026-10-07 — v2.0.0 — Outputs and input parameters optimized
#   2025-07-02 — v1.0.0 — Initial creation
#
###############################################################################


# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("--metadata", "-m"), type = "character", default = NULL,
              dest = "metadata",
              help = "Sample metadata file (TSV format)"),
  make_option(c("--expression", "-e"), type = "character", default = NULL,
              dest = "expression",
              help = "Count matrix file (TSV format)"),
  make_option(c("--gtf", "-g"), type = "character", default = NULL,
              dest = "gtf",
              help = "Annotation file (GTF) used to map transcript IDs to gene symbols"),
  make_option(c("--output-degs", "-d"), type = "character", default = "degs.tsv",
              dest = "output_degs",
              help = "Output file for differentially expressed genes [default %default]"),
  make_option(c("--output-plots", "-p"), type = "character",
              default = "pca_plot.png,heatmap.png",
              dest = "output_plots",
              help = "Output plot files, PCA and heatmap (comma-separated) [default %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Run DESeq2 differential expression analysis")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$metadata)) {
  stop("Error: --metadata argument is required")
}
if (is.null(opt$expression)) {
  stop("Error: --expression argument is required")
}
if (is.null(opt$gtf)) {
  stop("Error: --gtf argument is required")
}
if (!file.exists(opt$gtf)) {
  stop("Error: GTF file not found: ", opt$gtf)
}


# Read input files
cat("Reading input files...\n")
samples <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)
count_matrix <- read.table(opt$expression, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Ensure sample order matches between metadata and count matrix
samples <- samples[match(colnames(count_matrix), samples$sample), ]

# Create DESeq2 dataset
suppressWarnings({
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = samples,
                                design = ~ condition)
})

# Set reference level
dds$condition <- relevel(dds$condition, ref = "normal")

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Retained", nrow(dds), "genes for analysis\n")

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
suppressWarnings({
  suppressMessages({
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "cancer", "normal"))
  })
})

# Extract gene symbols from transcript/gene IDs
extract_gene_symbols <- function(transcript_ids, gtf_file) {
  gene_symbols <- rep(NA_character_, length(transcript_ids))

  cat("Extracting gene symbols from", gtf_file, "...\n")

  gtf_lines <- readLines(gtf_file)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

  # Create a lookup table for transcript_id to gene_name
  transcript_to_gene <- list()

  for (line in gtf_lines) {
    transcript_match <- regexpr('transcript_id "([^"]+)"', line)
    gene_name_match <- regexpr('gene_name "([^"]+)"', line)

    if (transcript_match > 0 && gene_name_match > 0) {
      transcript_id <- regmatches(line, transcript_match)
      transcript_id <- gsub('transcript_id "([^"]+)"', '\\1', transcript_id)

      gene_name <- regmatches(line, gene_name_match)
      gene_name <- gsub('gene_name "([^"]+)"', '\\1', gene_name)

      # Remove version number from transcript ID for matching
      clean_transcript_id <- sub("\\.[0-9]+$", "", transcript_id)
      transcript_to_gene[[clean_transcript_id]] <- gene_name
    }
  }

  # Map transcript IDs to gene symbols
  for (i in seq_along(transcript_ids)) {
    clean_id <- sub("\\.[0-9]+$", "", transcript_ids[i])

    if (clean_id %in% names(transcript_to_gene)) {
      gene_symbols[i] <- transcript_to_gene[[clean_id]]
    }
  }

  cat("Mapped", sum(!is.na(gene_symbols)), "out of", length(transcript_ids),
      "transcripts to gene symbols\n")

  # If no symbol found, use the transcript ID without version
  gene_symbols[is.na(gene_symbols)] <- sub("\\.[0-9]+$", "", transcript_ids[is.na(gene_symbols)])

  return(gene_symbols)
}


# Create output table with gene symbols
res_df <- data.frame(
  transcript_id = rownames(res),
  gene_symbol = extract_gene_symbols(rownames(res), opt$gtf),
  log2FoldChange = res$log2FoldChange,
  padj = res$padj,
  stringsAsFactors = FALSE
)

# Remove rows with NA padj values and sort by padj
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]

# Count up- and downregulated genes (padj < 0.05, |log2FC| > 1)
n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1)
cat("DEGs (padj < 0.05, |log2FC| > 1):", n_up + n_down, "\n")
cat("  Upregulated in cancer:  ", n_up, "\n")
cat("  Downregulated in cancer:", n_down, "\n")

# Write DEGs table
write.table(res_df, file = opt$output_degs, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Differentially expressed genes written to:", opt$output_degs, "\n")

# Parse plot output filenames
plot_files <- strsplit(opt$output_plots, ",")[[1]]
if (length(plot_files) != 2) {
  stop("Error: Exactly 2 plot files must be specified (PCA, heatmap)")
}

pca_file <- plot_files[1]
heatmap_file <- plot_files[2]

# Generate plots
cat("Generating plots...\n")

# Variance stabilizing transformation
suppressWarnings({
  suppressMessages({
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  })
})

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

png(pca_file, width = 800, height = 600)
plot(pca_data$PC1, pca_data$PC2,
     main = "PCA: Cancer vs Normal Tissue",
     xlab = paste0("PC1: ", percent_var[1], "% variance"),
     ylab = paste0("PC2: ", percent_var[2], "% variance"),
     col = ifelse(pca_data$condition == "cancer", "red", "blue"),
     pch = 16, cex = 2)
text(pca_data$PC1, pca_data$PC2, labels = pca_data$name, pos = 3, cex = 0.8)
legend("topright", legend = c("Cancer", "Normal"), col = c("red", "blue"), pch = 16)
dev.off()

# Heatmap of the top 20 genes by adjusted p-value, labelled with gene symbols
top_ids <- head(res_df$transcript_id, 20)
top_ids <- top_ids[top_ids %in% rownames(vsd)]
cat("Selected", length(top_ids), "genes for heatmap\n")

if (length(top_ids) > 1) {
  heatmap_data <- assay(vsd)[top_ids, , drop = FALSE]
  rownames(heatmap_data) <- make.unique(
    res_df$gene_symbol[match(top_ids, res_df$transcript_id)]
  )

  # Row-wise z-scores
  heatmap_z <- t(scale(t(heatmap_data)))
  heatmap_z[is.na(heatmap_z)] <- 0

  # Hierarchical clustering of rows and columns (reversed rows: first cluster on top)
  row_ord <- rev(hclust(dist(heatmap_z))$order)
  col_ord <- hclust(dist(t(heatmap_z)))$order

  heatmap_z <- heatmap_z[row_ord, col_ord, drop = FALSE]
  heat_cols <- colorRampPalette(c("blue", "white", "red"))(50)
  zmax <- max(abs(heatmap_z))

  png(heatmap_file, width = 900, height = 700)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))

  # Heatmap panel
  par(mar = c(8, 9, 3, 1))
  image(x = seq_len(ncol(heatmap_z)), y = seq_len(nrow(heatmap_z)),
        z = t(heatmap_z), col = heat_cols, zlim = c(-zmax, zmax),
        axes = FALSE, xlab = "", ylab = "",
        main = paste("Top", nrow(heatmap_z), "genes by adjusted p-value"))
  axis(1, at = seq_len(ncol(heatmap_z)), labels = colnames(heatmap_z), las = 2)
  axis(2, at = seq_len(nrow(heatmap_z)), labels = rownames(heatmap_z), las = 1, cex.axis = 0.8)
  box()

  # Color legend panel
  par(mar = c(8, 1, 3, 5))
  legend_z <- seq(-zmax, zmax, length.out = 50)
  image(x = 1, y = legend_z, z = matrix(legend_z, nrow = 1), col = heat_cols,
        axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1)
  mtext("Row z-score (VST counts)", side = 4, line = 3)
  box()

  dev.off()
  cat("Heatmap saved to:", heatmap_file, "\n")
} else {
  cat("Error: Fewer than 2 genes available for heatmap\n")
}

cat("Analysis complete. Plots saved as:", paste(plot_files, collapse = ", "), "\n")

