#!/bin/bash

# Input validation
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <star_bam_file> <output_directory>"
    echo "Example: $0 /path/to/Aligned.sortedByCoord.out.bam /path/to/output_dir"
    exit 1
fi

# Input parameters
BAM_FILE=$1
OUTPUT_DIR=$2

# Validate input file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Paths to required files - adjust these as needed
ARRIBA_DIR="/home/viagen/rnaseq/arriba_v2.4.0"
GENOME_FASTA="/home/viagen/rnaseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/home/viagen/rnaseq/Homo_sapiens.GRCh38.102.chr.gtf"

zcat "$ARRIBA_DIR/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz" > combined_fusions.tsv
cat "/home/viagen/rnaseq/arriba_v2.4.0/database/known_fusions_Gokce.tsv" >> combined_fusions.tsv

# Run Arriba without blacklist
echo "Running analysis without blacklist..."
  $ARRIBA_DIR/arriba \
    -x "$BAM_FILE" \
    -o "${OUTPUT_DIR}/fusions_all_events3.tsv" \
    -a "$GENOME_FASTA" \
    -g "$GTF_FILE" \
    -k combined_fusions.tsv \
    -t $ARRIBA_DIR/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
    -p $ARRIBA_DIR/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
    -S 1 \
    -E 5 \
    -l 200 \
    -A 15 \
    -m 0.9 \
    -L 0.4 \
    -R 150000 \
    -e 0.9 \
    -V 0.02 \
    -F 500 \
    -C 0.02 \
    -f same_gene,blacklist

echo "Analysis complete! Output file: ${OUTPUT_DIR}/fusions_noblacklist.tsv"