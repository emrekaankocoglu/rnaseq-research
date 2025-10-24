#!/bin/bash

export STAR_DIR="/home/viagen/rnaseq/STAR-2.7.11b/bin/Linux_x86_64_static"
export PATH="$STAR_DIR:$PATH"

GENOME_FASTA="/home/viagen/rnaseq/GRCh38.primary_assembly.genome.fa" #gencode bu 
GTF_FILE="/home/viagen/rnaseq/gencode.v44.basic.annotation.gtf"
GENOME_DIR="/home/viagen/rnaseq/star_genome_gencode"

$STAR_DIR/STAR --runMode genomeGenerate \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $GTF_FILE \
    --runThreadN 8 \

    