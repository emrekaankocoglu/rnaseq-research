#!/usr/bin/env Rscript

# Load required libraries
library(DT)
library(dplyr)
library(magrittr)
library(htmlwidgets)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript create_fusion_table.R <input_tsv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

# Extract patient ID and fusion type from filename
filename <- basename(input_file)
fusion_type <- ifelse(grepl("blacklist", filename), "blacklisted", "non-blacklisted")
ref_type <- ifelse(grepl("ensembl", filename), "Ensembl", "Gencode")

# Check if file exists
if (!file.exists(input_file)) {
  stop(paste("No fusion data found:", input_file))
}

# Read the data
fusions <- read.delim(input_file)

# Clean confidence levels for better sorting
confidence_order <- c("high", "medium", "low")
fusions$confidence <- factor(fusions$confidence, levels = confidence_order)

# Create and save the table
fusion_table <- datatable(fusions,
  options = list(
    pageLength = 25,
    scrollX = TRUE,
    dom = '<"top"f>rt<"bottom"lip>',
    order = list(list(13, 'desc')),
    columnDefs = list(
      list(
        targets = "_all",
        render = JS(
          "function(data, type, row, meta) {",
          "  return data === null ? '-' : data;",
          "}"
        )
      )
    )
  ),
  class = 'cell-border stripe',
  filter = 'top',
  rownames = FALSE,
  caption = htmltools::tags$caption(
    style = 'caption-side: top; text-align: center; font-size: 16px; font-weight: bold;',
    sprintf('%s Fusion Analysis (%s)', ref_type, fusion_type)
  )
) %>%
formatStyle(
  columns = 'confidence',
  backgroundColor = styleEqual(
    confidence_order,
    c('#FFFFFF', '#FFFFFF', '#FFFFFF')
  ),
  fontWeight = styleEqual(
    confidence_order,
    c('bold', 'normal', 'normal')
  )
)

# Create output filename
output_file <- file.path(output_dir, sub("\\.tsv$", "_table.html", basename(input_file)))

# Save as HTML file
saveWidget(fusion_table, output_file, selfcontained = TRUE)