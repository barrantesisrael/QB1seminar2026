#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("--input", "-i"), type="character", default="salmon_quant/",
              help="Input directory containing salmon quantification results [default %default]"),
  make_option(c("--metadata", "-m"), type="character", default="samples.tsv",
              help="Sample metadata file (TSV format) [default %default]"),
  make_option(c("--output", "-o"), type="character", default="counts.tsv",
              help="Output file for count matrix [default %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="Load salmon quantification data and create count matrix")
opt <- parse_args(opt_parser)

# Read sample metadata
if (!file.exists(opt$metadata)) {
  stop("Error: metadata file ", opt$metadata, " not found")
}

samples <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Ensure we only use the sample column
sample_names <- as.character(samples$sample)

# Construct file paths for salmon quantification files
files <- file.path(opt$input, sample_names, "quant.sf")
names(files) <- sample_names

# Check if all quantification files exist
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Error: Missing quantification files: ", paste(missing_files, collapse=", "))
}

# Read salmon quantification data and create count matrix
cat("Reading salmon quantification files...\n")

# Read first file to get gene/transcript names
first_data <- read.table(files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
gene_names <- first_data$Name

# Initialize count matrix
count_matrix <- matrix(0, nrow=length(gene_names), ncol=length(files))
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- sample_names

# Read all quantification files
for (i in 1:length(files)) {
  data <- read.table(files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
  count_matrix[, i] <- round(data$NumReads)
}

# Write count matrix to output file
write.table(count_matrix, file=opt$output, sep="\t", quote=FALSE, 
            row.names=TRUE, col.names=NA)

cat("Count matrix written to:", opt$output, "\n")
cat("Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")