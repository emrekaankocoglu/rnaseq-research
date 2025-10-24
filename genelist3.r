# Load required libraries
library(pathview)
library(org.Hs.eg.db)
library(KEGGREST)
library(clusterProfiler)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript genelist3.r <significant_genes_file> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

# WITH THIS:
pathway_dir <- file.path(output_dir, "expression_results", "pathways")
dir.create(pathway_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Read the TSV file
cat("Reading expression data...\n")
expr_data <- read.delim(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Step 2: Filter genes with Log2FC between abs(2) and abs(5)
filtered_data <- expr_data[abs(expr_data$Log2FC) >= 2 &
                          abs(expr_data$Log2FC) <= 5, ]

cat("Found", nrow(filtered_data), "genes with Log2FC between 2 and 5\n")

# Step 3: Get ENTREZ IDs for the filtered genes
entrez_ids <- mapIds(org.Hs.eg.db,
                    keys = filtered_data$external_gene_name,
                    keytype = "SYMBOL",
                    column = "ENTREZID")

# Remove any NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

if(length(entrez_ids) == 0) {
    cat("No valid gene symbols found. Exiting...\n")
    quit()
}

# Step 4: Find KEGG pathways
cat("\nSearching for KEGG pathways...\n")
kegg_pathways <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05
)

if(length(kegg_pathways@result$ID) == 0) {
    cat("No significant pathways found.\n")
    quit()
}

# Step 5: Create pathway visualizations
cat("\nGenerating pathway visualizations...\n")

# Create named vector of fold changes
gene_fc <- filtered_data$Log2FC
names(gene_fc) <- mapIds(org.Hs.eg.db,
                        keys = filtered_data$external_gene_name,
                        keytype = "SYMBOL",
                        column = "ENTREZID")

# Remove any NA values
gene_fc <- gene_fc[!is.na(names(gene_fc))]

setwd(output_dir)
# REPLACE THE pathview LOOP WITH:
for(pathway in kegg_pathways@result$ID) {
    tryCatch({
        pathview(
            gene.data = gene_fc,
            pathway.id = pathway,
            species = "hsa",
            kegg.native = TRUE,
            limit = list(gene = max(abs(gene_fc))),
            low = list(gene = "blue"),
            mid = list(gene = "white"),
            high = list(gene = "red"),
            out.suffix = paste0("pathway_", pathway),
            filename = paste0("pathway_", pathway)
        )
    }, error = function(e) {
        cat("Error processing pathway:", pathway, "\n")
        cat("Error message:", e$message, "\n")
    })
}

write.csv(as.data.frame(kegg_pathways),
          file.path(pathway_dir, "pathway_enrichment_results.csv"),
          row.names = FALSE)
write.csv(filtered_data,
          file.path(pathway_dir, "filtered_genes.csv"),
          row.names = FALSE)


# Print summary
cat("\nAnalysis complete! Results saved in", pathway_dir, "folder:\n")
cat("- KEGG pathway visualizations\n")
cat("- Enrichment results in 'pathway_enrichment_results.csv'\n")
cat("- Filtered genes in 'filtered_genes.csv'\n")

# Print top pathways
cat("\nTop enriched pathways:\n")
print(head(kegg_pathways))