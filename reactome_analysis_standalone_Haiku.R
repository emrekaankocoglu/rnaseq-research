#!/usr/bin/env Rscript

# Command line argument handling
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript reactome_analysis_standalone.R <tcga_reference_file> <star_output_file> <patient_name>")
}

tcga_reference_file <- args[1]
star_output_file <- args[2]
patient_name <- args[3]

# Load required libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(biomaRt)
    library(knitr)
    library(kableExtra)
})

# Define paths
output_dir <- dirname(tcga_reference_file)
results_dir <- file.path(output_dir, "expression_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read TCGA reference data
message("Reading TCGA reference data...")
tcga_data <- readLines(tcga_reference_file)
tcga_sample_id <- tcga_data[1]
tcga_counts <- data.frame(do.call(rbind, strsplit(tcga_data[-1], "\t")))
colnames(tcga_counts) <- c("gene_id", "count")
tcga_counts$count <- as.numeric(tcga_counts$count)
tcga_counts$gene_id <- sub("\\..*", "", tcga_counts$gene_id)

# Read STAR output
message("Reading STAR output...")
star_counts <- read.table(star_output_file, skip=4)
colnames(star_counts) <- c("gene_id", "unstranded", "first_strand", "second_strand")
your_counts <- data.frame(
    gene_id = star_counts$gene_id,
    counts = star_counts$second_strand
)
your_counts$gene_id <- sub("\\..*", "", your_counts$gene_id)

# Find common genes
common_genes <- intersect(your_counts$gene_id, tcga_counts$gene_id)
message("Number of common genes found: ", length(common_genes))

# Create comparison table
your_data_subset <- your_counts[match(common_genes, your_counts$gene_id), ]
tcga_data_subset <- tcga_counts[match(common_genes, tcga_counts$gene_id), ]

comparison <- data.frame(
    Gene = your_data_subset$gene_id,
    Reference_Count = tcga_data_subset$count,
    Your_Count = your_data_subset$counts,
    stringsAsFactors = FALSE
)

# Calculate fold changes
comparison$Log2FC <- log2((comparison$Your_Count + 1) / (comparison$Reference_Count + 1))

# Add gene names
message("Adding gene annotations...")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id",
    values = comparison$Gene,
    mart = ensembl
)

# Merge gene info with comparison data
comparison <- merge(comparison, gene_info, 
                   by.x = "Gene", 
                   by.y = "ensembl_gene_id", 
                   all.x = TRUE)

# Define significance thresholds
fc_threshold <- 1.5

# Add regulation status
comparison$Status <- "No Change"
comparison$Status[comparison$Log2FC > log2(fc_threshold)] <- "Up"
comparison$Status[comparison$Log2FC < -log2(fc_threshold)] <- "Down"

# Sort by absolute fold change
comparison <- comparison[order(-abs(comparison$Log2FC)),]

# Write complete comparison table
write.table(comparison, 
            file=file.path(results_dir, "complete_expression_analysis.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create significant genes summary
sig_genes <- comparison[comparison$Status != "No Change",]
sig_genes <- sig_genes[order(-abs(sig_genes$Log2FC)),]

# Write significant genes table
write.table(sig_genes,
            file=file.path(results_dir, "significant_genes.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Filter for 3-fold changes (log2(3) ≈ 1.58)
fold_3_changes <- comparison[abs(comparison$Log2FC) > 1.58, ]
fold_3_changes <- fold_3_changes[order(-abs(fold_3_changes$Log2FC)), ]

# Separate up and down regulated for 3-fold changes
up_regulated_3fold <- fold_3_changes[fold_3_changes$Log2FC > 1.58, ]
down_regulated_3fold <- fold_3_changes[fold_3_changes$Log2FC < -1.58, ]

# Sort both by absolute fold change
up_regulated_3fold <- up_regulated_3fold[order(-up_regulated_3fold$Log2FC), ]
down_regulated_3fold <- down_regulated_3fold[order(down_regulated_3fold$Log2FC), ]

# Write 3-fold change files
write.table(up_regulated_3fold[, c("Gene", "external_gene_name", "Log2FC")],
            file=file.path(results_dir, "threefold_up.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

write.table(down_regulated_3fold[, c("Gene", "external_gene_name", "Log2FC")],
            file=file.path(results_dir, "threefold_down.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create Reactome input files
reactome_input <- data.frame(
    Gene = comparison$Gene,
    FoldChange = comparison$Log2FC
)

write.table(reactome_input,
            file=file.path(results_dir, "reactome_input.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create separate files for Reactome analysis
write.table(sig_genes[sig_genes$Status == "Up", c("Gene", "external_gene_name", "Log2FC")],
            file=file.path(results_dir, "upregulated_genes.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

write.table(sig_genes[sig_genes$Status == "Down", c("Gene", "external_gene_name", "Log2FC")],
            file=file.path(results_dir, "downregulated_genes.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create comprehensive summary report
sink(file.path(results_dir, "analysis_summary.txt"))
cat("=== Analysis Summary for", patient_name, "===\n\n")
cat("TCGA Sample ID:", tcga_sample_id, "\n")
cat("Total genes analyzed:", nrow(comparison), "\n")
cat("\n=== Significant Changes (1.5-fold) ===\n")
cat("Significantly upregulated genes:", sum(sig_genes$Status == "Up"), "\n")
cat("Significantly downregulated genes:", sum(sig_genes$Status == "Down"), "\n")

cat("\n=== 3-Fold Changes ===\n")
cat("Total genes with >3-fold change:", nrow(fold_3_changes), "\n")
cat("Up-regulated (>3-fold):", nrow(up_regulated_3fold), "\n")
cat("Down-regulated (<-3-fold):", nrow(down_regulated_3fold), "\n")

cat("\nTop 10 Upregulated Genes (>3-fold):\n")
print(head(up_regulated_3fold[, c("external_gene_name", "Log2FC")], 10))

cat("\nTop 10 Downregulated Genes (<-3-fold):\n")
print(head(down_regulated_3fold[, c("external_gene_name", "Log2FC")], 10))
sink()

# Print final status message
message("\nAnalysis complete! Files created in: ", results_dir)
message("\nKey output files:")
message("1. Complete analysis: ", file.path(results_dir, "complete_expression_analysis.tsv"))
message("2. Significant genes: ", file.path(results_dir, "significant_genes.tsv"))
message("3. Three-fold up regulated: ", file.path(results_dir, "threefold_up.tsv"))
message("4. Three-fold down regulated: ", file.path(results_dir, "threefold_down.tsv"))
message("5. Reactome input: ", file.path(results_dir, "reactome_input.txt"))
message("6. Analysis summary: ", file.path(results_dir, "analysis_summary.txt"))