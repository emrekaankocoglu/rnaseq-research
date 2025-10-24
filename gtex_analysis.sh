#!/bin/bash

PATIENT_NAME=$1
READSPERGENE_FILE=$2
OUTPUT_DIR=$3
REF_TYPE=$4
VERSION_INFO=$5
READSPERGENE_FILENAME=$6

# Create R script for GTEx analysis
cat > "$OUTPUT_DIR/gtex_select.R" << 'EOF'
suppressPackageStartupMessages({
    library(recount3)
    library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
mode <- args[2]
tissue_num <- if(length(args) >= 3) as.numeric(args[3]) else NULL
sample_num <- if(length(args) >= 4) as.numeric(args[4]) else NULL

# Get tissues
projects <- available_projects()
gtex_tissues <- subset(projects, project_type == "data_sources" & file_source == "gtex")
gtex_tissues <- gtex_tissues[order(gtex_tissues$project), ]

if (mode == "list_tissues") {
    for(i in seq_len(nrow(gtex_tissues))) {
        tissue_name <- gsub("^gtex_", "", gtex_tissues$project[i])
        cat(sprintf("%d. %s (%d samples)\n", i, tissue_name, gtex_tissues$n_samples[i]))
    }
    saveRDS(gtex_tissues, file=file.path(output_dir, "gtex_tissues.rds"))

} else if (mode == "list_samples") {
    gtex_tissues <- readRDS(file.path(output_dir, "gtex_tissues.rds"))
    selected_tissue <- gtex_tissues[tissue_num, ]
    
    tryCatch({
        rse_gene <- create_rse(selected_tissue)
        metadata <- as.data.frame(colData(rse_gene))
        
        for(i in seq_len(ncol(rse_gene))) {
            age <- if(is.na(metadata$gtex.age[i])) "Age:NA" else paste0("Age:", metadata$gtex.age[i])
            sex <- if(is.na(metadata$gtex.sex[i])) "Sex:NA" else paste0("Sex:", metadata$gtex.sex[i])
            site <- if(is.na(metadata$gtex.smtsd[i])) "Site:NA" else paste0("Site:", metadata$gtex.smtsd[i])
            quality <- if(is.na(metadata$gtex.smtsisch[i])) "QC:NA" else paste0("QC:", metadata$gtex.smtsisch[i])
            
            cat(sprintf("%d. %s | %s | %s | %s | %s\n",
                i, rownames(metadata)[i], age, sex, site, quality))
        }
        
        saveRDS(rse_gene, file=file.path(output_dir, "selected_tissue_data.rds"))
        
    }, error = function(e) {
        cat(sprintf("\nError accessing tissue data: %s\n", e$message))
        quit(status = 1)
    })

} else if (mode == "get_sample") {
    rse_gene <- readRDS(file.path(output_dir, "selected_tissue_data.rds"))
    gtex_tissues <- readRDS(file.path(output_dir, "gtex_tissues.rds"))
    
    metadata <- as.data.frame(colData(rse_gene))
    sample_id <- rownames(metadata)[sample_num]
    
    counts <- assay(rse_gene[, sample_id], "raw_counts")
    rownames(counts) <- rowData(rse_gene)$gene_id
    
    write.table(
        data.frame(gene_id=rownames(counts), count=counts),
        file=file.path(output_dir, "normal_reference_counts.tab"),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
    
    selected_info <- list(
        tissue = gtex_tissues$project[tissue_num],
        sample_id = sample_id,
        age = metadata$gtex.age[sample_num],
        sex = metadata$gtex.sex[sample_num],
        body_site = metadata$gtex.smtsd[sample_num],
        tissue_quality = metadata$gtex.smtsisch[sample_num]
    )
    saveRDS(selected_info, file=file.path(output_dir, "selected_tissue.rds"))
}
EOF

# Create log directory
mkdir -p "$OUTPUT_DIR/logs"

# Step 1: List available tissues
echo "Fetching available GTEx tissues..."
Rscript "$OUTPUT_DIR/gtex_select.R" "$OUTPUT_DIR" "list_tissues"


if [ $? -ne 0 ]; then
    echo "Error: Failed to fetch tissue list"
    exit 1
fi

# Get tissue selection
read -p "Enter tissue number: " TISSUE_NUM

# Validate tissue number
if ! [[ "$TISSUE_NUM" =~ ^[0-9]+$ ]]; then
    echo "Error: Please enter a valid number"
    exit 1
fi

# Step 2: List samples for selected tissue
echo "Fetching samples for selected tissue..."
Rscript "$OUTPUT_DIR/gtex_select.R" "$OUTPUT_DIR" "list_samples" "$TISSUE_NUM"

if [ $? -ne 0 ]; then
    echo "Error: Failed to fetch sample list"
    exit 1
fi

# Get sample selection
read -p "Enter sample number: " SAMPLE_NUM

# Validate sample number
if ! [[ "$SAMPLE_NUM" =~ ^[0-9]+$ ]]; then
    echo "Error: Please enter a valid number"
    exit 1
fi

# Step 3: Get selected sample data
echo "Processing selected sample..."
Rscript "$OUTPUT_DIR/gtex_select.R" "$OUTPUT_DIR" "get_sample" "$TISSUE_NUM" "$SAMPLE_NUM"

if [ $? -ne 0 ]; then
    echo "Error: Failed to process selected sample"
    exit 1
fi

# Create reference info file
selected_info=$(Rscript -e "info <- readRDS('$OUTPUT_DIR/selected_tissue.rds'); cat(sprintf('Tissue: %s\nSample ID: %s\nAge: %s\nSex: %s\nBody Site: %s\nTissue Quality: %s', info\$tissue, info\$sample_id, info\$age, info\$sex, info\$body_site, info\$tissue_quality))")

cat > "$OUTPUT_DIR/reference_info.txt" << EOF
=== GTEx Reference Information ===
Date: $(date)
$selected_info
EOF

# Update version info
VERSION_INFO+="Expression Analysis - ${REF_TYPE} - GTEx Reference ($(Rscript -e "cat(readRDS('$OUTPUT_DIR/selected_tissue.rds')\$tissue)")) - Input: ${READSPERGENE_FILENAME}"
echo "$VERSION_INFO"