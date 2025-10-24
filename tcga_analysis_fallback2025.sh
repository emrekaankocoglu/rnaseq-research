#!/bin/bash

# Input parameters
PATIENT_NAME=$1
READSPERGENE_FILE=$2
OUTPUT_DIR=$3
REF_TYPE=$4
VERSION_INFO=$5
READSPERGENE_FILENAME=$6

echo "=== TCGA Reference Selection ==="

# First ask for TCGA cancer code
echo "Please enter the TCGA cancer code (e.g., BRCA, LUAD, etc.):"
read CANCER_TYPE

# Get available options using clinical data - with error handling for data mismatches
echo "Retrieving available options for ${CANCER_TYPE}..."
Rscript -e "
library(TCGAbiolinks)
library(methods)
cancer_code <- commandArgs(trailingOnly=TRUE)[1]

# Use tryCatch to handle potential errors
tryCatch({
    # First query to make sure the cancer code exists
    query <- GDCquery(
        project = cancer_code,
        data.category = 'Transcriptome Profiling', 
        data.type = 'Gene Expression Quantification',
        workflow.type = 'STAR - Counts'
    )
    
    # Get sample types from the query results
    samples <- getResults(query)
    sample_types <- unique(samples\$sample_type)
    
    # Get clinical data with safe error handling
    clinical <- tryCatch({
        clin_data <- GDCquery_clinic(cancer_code)
        # If data.frame dimensions don't match, simplify and ensure consistency
        if(is.data.frame(clin_data)) {
            # Keep only essential columns to avoid mismatches
            essential_cols <- c('submitter_id', 'tissue_or_organ_of_origin', 'primary_diagnosis', 'sample_type')
            existing_cols <- intersect(essential_cols, colnames(clin_data))
            clin_data <- clin_data[, existing_cols, drop=FALSE]
        }
        clin_data
    }, error = function(e) {
        # Create minimal clinical data if query fails
        message('Note: Using simplified clinical data due to API changes')
        data.frame(
            submitter_id = unique(samples\$cases),
            tissue_or_organ_of_origin = rep('Unknown', length(unique(samples\$cases))),
            primary_diagnosis = rep('Unknown', length(unique(samples\$cases))),
            stringsAsFactors = FALSE
        )
    })
    
    cat('\nAvailable Sample Types:\n')
    for(i in seq_along(sample_types)) {
        cat(i, ': ', sample_types[i], '\n')
    }
    
    # Handle potentially missing columns gracefully
    if('tissue_or_organ_of_origin' %in% colnames(clinical)) {
        tissue_types <- unique(clinical\$tissue_or_organ_of_origin)
        cat('\nAvailable Tissue Types:\n')
        for(i in seq_along(tissue_types)) {
            cat(i, ': ', tissue_types[i], '\n')
        }
    } else {
        tissue_types <- c('Unknown')
        cat('\nTissue type information unavailable\n')
    }
    
    if('primary_diagnosis' %in% colnames(clinical)) {
        diagnoses <- unique(clinical\$primary_diagnosis)
        cat('\nAvailable Primary Diagnoses:\n')
        for(i in seq_along(diagnoses)) {
            cat(i, ': ', diagnoses[i], '\n')
        }
    } else {
        diagnoses <- c('Unknown')
        cat('\nDiagnosis information unavailable\n')
    }
    
    # Save the options for later use
    selections <- list(
        sample_types = sample_types,
        tissue_types = tissue_types,
        diagnoses = diagnoses
    )
    saveRDS(selections, file=paste0(commandArgs(trailingOnly=TRUE)[2], '/tcga_selections.rds'))
}, error = function(e) {
    cat('\nError retrieving data:', conditionMessage(e), '\n')
    # Provide informative error to help debugging
    if(grepl('was not found', conditionMessage(e))) {
        cat('\nAvailable TCGA projects:\n')
        projects <- getGDCprojects()
        print(projects[,c('id', 'name')])
    }
    quit(status = 1)
})
" "$CANCER_TYPE" "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo "Error: Failed to retrieve TCGA data options. Please check the cancer code and try again."
    exit 1
fi

# Get user selections
read -p "Enter sample type number (or 'all' to see all cases): " sample_type_num
read -p "Enter tissue type number (or 'all' to see all cases): " tissue_type_num 
read -p "Enter primary diagnosis number (or 'all' to see all cases): " diagnosis_num

# Get matching samples with enhanced error handling
echo "Retrieving matching samples..."
Rscript -e "
library(TCGAbiolinks)
library(methods)
args <- commandArgs(trailingOnly=TRUE)
cancer_code <- args[1]
output_dir <- args[2]
sample_type_input <- args[3]
tissue_type_input <- args[4]
diagnosis_input <- args[5]

# Load saved selections
selections <- readRDS(paste0(output_dir, '/tcga_selections.rds'))

# Get RNA-seq data
tryCatch({
    query <- GDCquery(
        project = cancer_code,
        data.category = 'Transcriptome Profiling',
        data.type = 'Gene Expression Quantification', 
        workflow.type = 'STAR - Counts'
    )
    samples <- getResults(query)
    
    # Handle case filtering with fall-backs for missing data
    if(sample_type_input == 'all' || tissue_type_input == 'all' || diagnosis_input == 'all') {
        filtered_samples <- samples
    } else {
        # Convert inputs to numeric and handle possible NAs
        sample_type_idx <- as.numeric(sample_type_input)
        tissue_type_idx <- as.numeric(tissue_type_input)
        diagnosis_idx <- as.numeric(diagnosis_input)
        
        # Get selected values
        selected_sample_type <- selections\$sample_types[sample_type_idx]
        
        # Filter primarily by sample type which should always be available
        filtered_samples <- samples[samples\$sample_type == selected_sample_type,]
        
        # Apply additional filters only if data is available and not 'Unknown'
        if(tissue_type_input != 'all' && length(selections\$tissue_types) > 1 && 
           any(selections\$tissue_types != 'Unknown')) {
            selected_tissue <- selections\$tissue_types[tissue_type_idx]
            
            # Get clinical data with error handling
            clinical <- tryCatch({
                GDCquery_clinic(cancer_code)
            }, error = function(e) {
                message('Note: Using simplified data due to API changes')
                data.frame(
                    submitter_id = unique(samples\$cases),
                    stringsAsFactors = FALSE
                )
            })
            
            # Filter by tissue type if the column exists
            if('tissue_or_organ_of_origin' %in% colnames(clinical)) {
                matching_cases <- clinical\$submitter_id[
                    clinical\$tissue_or_organ_of_origin == selected_tissue
                ]
                filtered_samples <- filtered_samples[
                    filtered_samples\$cases %in% matching_cases,
                ]
            }
            
            # Filter by diagnosis if the column exists and diagnosis filter selected
            if(diagnosis_input != 'all' && 'primary_diagnosis' %in% colnames(clinical) &&
               length(selections\$diagnoses) > 1 && any(selections\$diagnoses != 'Unknown')) {
                selected_diagnosis <- selections\$diagnoses[diagnosis_idx]
                matching_cases <- clinical\$submitter_id[
                    clinical\$primary_diagnosis == selected_diagnosis
                ]
                filtered_samples <- filtered_samples[
                    filtered_samples\$cases %in% matching_cases,
                ]
            }
        }
    }
    
    if(nrow(filtered_samples) > 0) {
        cat('\nMatching Sample IDs:\n')
        for(i in seq_len(nrow(filtered_samples))) {
            cat(sprintf('%d. %s (%s)\n', i, 
                        filtered_samples\$cases[i], 
                        filtered_samples\$sample_type[i]))
        }
        saveRDS(filtered_samples, file=paste0(output_dir, '/filtered_samples.rds'))
    } else {
        cat('\nNo samples found matching these criteria.\n')
        quit(status = 1)
    }
}, error = function(e) {
    cat('\nError retrieving samples:', conditionMessage(e), '\n')
    quit(status = 1)
})
" "$CANCER_TYPE" "$OUTPUT_DIR" "$sample_type_num" "$tissue_type_num" "$diagnosis_num"

if [ $? -ne 0 ]; then
    echo "Error: Failed to process samples. Please try again with different criteria."
    exit 1
fi

# Get sample selection from user
read -p "Enter the number of the sample you want to use: " sample_num

# Download and prepare reference data
echo "Downloading and processing selected sample..."
Rscript -e "
library(TCGAbiolinks)
library(SummarizedExperiment)
library(methods)
args <- commandArgs(trailingOnly=TRUE)
output_dir <- args[1]
sample_idx <- as.numeric(args[2])
cancer_code <- args[3]

tryCatch({
    filtered_samples <- readRDS(paste0(output_dir, '/filtered_samples.rds'))
    selected_case <- na.omit(filtered_samples\$cases)[sample_idx]
    selected_sample_type <- na.omit(filtered_samples\$sample_type)[sample_idx]
    
    query <- GDCquery(
        project = cancer_code,
        data.category = 'Transcriptome Profiling',
        data.type = 'Gene Expression Quantification',
        workflow.type = 'STAR - Counts',
        sample.type = selected_sample_type,
        barcode = selected_case
    )
    
    GDCdownload(query, method = 'api', directory = paste0(output_dir, '/GDCdownload'))
    data <- GDCprepare(query, directory = paste0(output_dir, '/GDCdownload'))
    counts <- assay(data)
    
    write.table(counts, 
               file = paste0(output_dir, '/normal_reference_counts.tab'),
               sep = '\t', 
               quote = FALSE)
    
    write.table(
        data.frame(
            CANCER_TYPE = cancer_code,
            SAMPLE_ID = selected_case,
            SAMPLE_TYPE = selected_sample_type
        ),
        file = paste0(output_dir, '/sample_context.txt'),
        sep = '\t',
        quote = FALSE,
        row.names = FALSE
    )
}, error = function(e) {
    cat('\nError processing selected sample:', conditionMessage(e), '\n')
    quit(status = 1)
})
" "$OUTPUT_DIR" "$sample_num" "$CANCER_TYPE"

if [ $? -ne 0 ]; then
    echo "Error: Failed to download and process the selected sample."
    exit 1
fi

# Read sample context for VERSION_INFO
CONTEXT=$(cat "$OUTPUT_DIR/sample_context.txt" | tail -n 1)
VERSION_INFO+="Expression Analysis - ${REF_TYPE} - Reference ($CONTEXT)"
echo "$VERSION_INFO"
echo "Reference data preparation complete. Proceeding to differential expression analysis..."