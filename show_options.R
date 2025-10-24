library(TCGAbiolinks)

# Get arguments
args <- commandArgs(trailingOnly=TRUE)
cancer_code <- args[1]

# Get clinical data
clinical <- GDCquery_clinic(paste0("TCGA-", cancer_code))

# Display diagnosis types
cat("\nAvailable diagnosis types:\n")
diagnoses <- na.omit(unique(clinical$primary_diagnosis))
for(i in seq_along(diagnoses)) {
    cat(sprintf("%d. %s\n", i, diagnoses[i]))
}

# Display tissue types
cat("\nAvailable tissue types:\n")
tissues <- na.omit(unique(clinical$tissue_or_organ_of_origin))
for(i in seq_along(tissues)) {
    cat(sprintf("%d. %s\n", i, tissues[i]))
}

# Display sample types
cat("\nAvailable sample types:\n")
samples <- na.omit(unique(clinical$sample_type))
for(i in seq_along(samples)) {
    cat(sprintf("%d. %s\n", i, samples[i]))
}

# Save options
saveRDS(list(diagnoses=diagnoses, tissues=tissues, samples=samples), "options.rds")
