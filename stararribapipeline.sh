#!/bin/bash
#1st input validation: choice 1-2 
validate_input() {   
    local input=$1    
    local min=$2     
    local max=$3     

    if ! [[ "$input" =~ ^[0-9]+$ ]] ||    
       [ "$input" -lt "$min" ] ||         
       [ "$input" -gt "$max" ]; then      
        echo "Invalid selection. Please choose a number between $min and $max."
        exit 1
    fi
}

BASE_DIR="/home/viagen/rnaseq"
export STAR_DIR="/home/viagen/rnaseq/STAR-2.7.11b/bin/Linux_x86_64_static"
export ARRIBA_DIR="/home/viagen/rnaseq/arriba_v2.4.0"
export PATH="$STAR_DIR:$ARRIBA_DIR:$PATH"
export TRIMMOMATIC_DIR="${BASE_DIR}/Trimmomatic-0.39"
export BBMAP_DIR="${BASE_DIR}/bbmap"

ENSEMBL_DIR="/home/viagen/rnaseq/star_genome"
GENCODE_DIR="/home/viagen/rnaseq/star_genome_gencode"
ENSEMBL_FASTA="/home/viagen/rnaseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENCODE_FASTA="/home/viagen/rnaseq/GRCh38.primary_assembly.genome.fa"
ENSEMBL_GTF="/home/viagen/rnaseq/Homo_sapiens.GRCh38.102.chr.gtf"
GENCODE_GTF="/home/viagen/rnaseq/gencode.v44.basic.annotation.gtf"

log_and_update_history() {
    local message="$1"
    local log_file="$2"
    local history_file="$3"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$log_file"
    echo "$VERSION_INFO - $message" >> "$history_file"
}

# Add cleanup function and traps
cleanup() {
    local exit_code=$?
    local signal=$1
    
    if [ $exit_code -ne 0 ]; then
        log_and_update_history "Pipeline interrupted by signal $signal (exit code: $exit_code)" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
    fi
    
    sync
    
    exit $exit_code
}

# Set up traps
trap 'cleanup SIGHUP' SIGHUP
trap 'cleanup SIGINT' SIGINT
trap 'cleanup SIGTERM' SIGTERM

run_qc_pipeline() {
    local fastq_1=$1
    local fastq_2=$2
    local output_dir=$3
    
    echo "Starting QC filtering..."
    local qc_dir="${output_dir}/qc"
    mkdir -p "$qc_dir"

    # Step 1: Detect adapters using BBMerge with more reads sampled
    echo "Detecting adapters..."
    "${BBMAP_DIR}/bbmerge.sh" \
        in1="$fastq_1" \
        in2="$fastq_2" \
        outa="${qc_dir}/detected_adapters.fa" \
        reads=2m

    # Log the detected adapters
    echo "=== Detected Adapters ===" > "${qc_dir}/adapter_detection_log.txt"
    echo "Date: $(date)" >> "${qc_dir}/adapter_detection_log.txt"
    echo "Input files:" >> "${qc_dir}/adapter_detection_log.txt"
    echo "R1: $(basename "$fastq_1")" >> "${qc_dir}/adapter_detection_log.txt"
    echo "R2: $(basename "$fastq_2")" >> "${qc_dir}/adapter_detection_log.txt"
    echo -e "\nDetected adapter sequences:" >> "${qc_dir}/adapter_detection_log.txt"
    cat "${qc_dir}/detected_adapters.fa" >> "${qc_dir}/adapter_detection_log.txt"

    # Step 2: More stringent trimming with Trimmomatic
    echo "Running Trimmomatic..."
    java -jar "${TRIMMOMATIC_DIR}/trimmomatic-0.39.jar" PE \
        "$fastq_1" "$fastq_2" \
        "${qc_dir}/filtered_R1.fastq.gz" "${qc_dir}/unpaired_R1.fastq.gz" \
        "${qc_dir}/filtered_R2.fastq.gz" "${qc_dir}/unpaired_R2.fastq.gz" \
        ILLUMINACLIP:"${qc_dir}/detected_adapters.fa":2:20:7:1:true \
        SLIDINGWINDOW:4:20 \
        LEADING:20 \
        TRAILING:20 \
        MINLEN:36 \
        AVGQUAL:25

    # Step 3: Remove PhiX and filter by GC content in one BBDuk command
    echo "Running BBDuk contamination filtering..."
    "${BBMAP_DIR}/bbduk.sh" \
        in1="${qc_dir}/filtered_R1.fastq.gz" \
        in2="${qc_dir}/filtered_R2.fastq.gz" \
        out1="${qc_dir}/cleaned_R1.fastq.gz" \
        out2="${qc_dir}/cleaned_R2.fastq.gz" \
        ref="${BBMAP_DIR}/resources/phix174_ill.ref.fa.gz" \
        k=31 \
        hdist=1 \
        mingc=0.2 \
        maxgc=0.7 \
        stats="${qc_dir}/contaminants_stats.txt"

    # Generate QC report
    echo "=== QC Processing Report ===" > "${qc_dir}/qc_report.txt"
    echo "Date: $(date)" >> "${qc_dir}/qc_report.txt"
    echo -e "\nInput Files:" >> "${qc_dir}/qc_report.txt"
    echo "R1: $(basename "$fastq_1")" >> "${qc_dir}/qc_report.txt"
    echo "R2: $(basename "$fastq_2")" >> "${qc_dir}/qc_report.txt"
    echo -e "\nOutput Files:" >> "${qc_dir}/qc_report.txt"
    echo "Cleaned R1: cleaned_R1.fastq.gz" >> "${qc_dir}/qc_report.txt"
    echo "Cleaned R2: cleaned_R2.fastq.gz" >> "${qc_dir}/qc_report.txt"
    
    echo "${qc_dir}/cleaned_R1.fastq.gz ${qc_dir}/cleaned_R2.fastq.gz"
}

# Get patient name
read -p "Enter patient name: " PATIENT_NAME

# Clean up patient name for directory naming
PATIENT_NAME=$(echo "$PATIENT_NAME" | tr ' ' '_' | tr -cd '[:alnum:]_-')

# Set up versioning
PATIENT_BASE_DIR="/home/viagen/rnaseq/outputfolder/${PATIENT_NAME}"
mkdir -p "$PATIENT_BASE_DIR"

# Version control function
get_next_version() {
    local patient_dir="$1"
    
   # Always look for highest existing directory number and add 1
    local highest_ver=$(ls -d "$patient_dir"/v* 2>/dev/null | grep -o '[0-9]\+' | sort -n | tail -1)
    if [ -z "$highest_ver" ]; then
        echo "1"
    else
        echo "$((highest_ver + 1))"
    fi
}

# Create or update version history
if [ ! -f "$PATIENT_BASE_DIR/version_history.txt" ]; then
    cat > "$PATIENT_BASE_DIR/version_history.txt" << EOF
=== Analysis Version History for $PATIENT_NAME ===
Created: $(date)

Format:
v<VERSION> - <DATE> - <ANALYSIS_TYPE> - <REFERENCE> - <DETAILS>

EOF
fi

# Get next version and create directory
VERSION=$(get_next_version "$PATIENT_BASE_DIR")
OUTPUT_DIR="$PATIENT_BASE_DIR/v${VERSION}"
mkdir -p "$OUTPUT_DIR"

echo "What would you like to do?"
echo "1) Run new analysis"
echo "2) Visualize specific fusion(s) from existing results"
read -p "Enter your choice (1 or 2): " initial_choice
validate_input "$initial_choice" 1 2

if [ "$initial_choice" = "2" ]; then
  # Get patient name and find existing files
    read -p "Enter patient name: " PATIENT_NAME
    PATIENT_NAME=$(echo "$PATIENT_NAME" | tr ' ' '_' | tr -cd '[:alnum:]_-')
    
    # List all available versions and their reference types
    echo -e "\nAvailable analysis versions for $PATIENT_NAME:"
    for version_dir in "${BASE_DIR}/outputfolder/${PATIENT_NAME}"/v*; do
        if [ -d "$version_dir" ]; then
            version_num=$(basename "$version_dir" | grep -o '[0-9]\+')
            
            # Check for ensembl or gencode in fusion files
            if [ -f "$version_dir/ensembl_fusions_noblacklist.tsv" ]; then
                echo "Version $version_num - Ensembl analysis"
            elif [ -f "$version_dir/gencode_fusions_noblacklist.tsv" ]; then
                echo "Version $version_num - Gencode analysis"
            fi
        fi
    done

    # Ask user which version to use
    read -p "Enter the version number you want to analyze: " version_choice
    PATIENT_DIR="${BASE_DIR}/outputfolder/${PATIENT_NAME}/v${version_choice}"
    
    if [ ! -d "$PATIENT_DIR" ]; then
        echo "Error: Version $version_choice not found for patient ${PATIENT_NAME}"
        exit 1
    fi

    # Determine reference type from existing files
    if [ -f "$PATIENT_DIR/ensembl_fusions_noblacklist.tsv" ]; then
        REF_TYPE="ensembl"
        FUSION_FILE="$PATIENT_DIR/ensembl_fusions_noblacklist.tsv"
        BAM_FILE="$PATIENT_DIR/ensembl_fusionaligned.Aligned.sortedByCoord.out.bam"
        GTF_FILE="$ENSEMBL_GTF"
        echo "Using Ensembl analysis results (optimized for cancer fusions and ITDs)"
    elif [ -f "$PATIENT_DIR/gencode_fusions_noblacklist.tsv" ]; then
        REF_TYPE="gencode"
        FUSION_FILE="$PATIENT_DIR/gencode_fusions_noblacklist.tsv"
        BAM_FILE="$PATIENT_DIR/gencode_fusionaligned.Aligned.sortedByCoord.out.bam"
        GTF_FILE="$GENCODE_GTF"
        echo "Using Gencode analysis results (optimized for Ig/TcR fusions)"
    else
        echo "Error: No fusion analysis results found in version $version_choice"
        exit 1
    fi
    
    if [ ! -f "$FUSION_FILE" ] || [ ! -f "$BAM_FILE" ]; then
        echo "Error: Required fusion analysis files not found in version $version_choice"
        exit 1
    fi

while true; do
    echo -e "\nPlease specify the fusion you want to visualize:"
    read -p "Enter gene1 name: " gene1
    read -p "Enter gene2 name: " gene2
    read -p "Enter first direction (e.g., +/+ or +/-): " direction1
    read -p "Enter second direction (e.g., +/+ or +/-): " direction2
    read -p "Enter confidence level (high/medium/low): " confidence

    # Create unique filename for this fusion
    FUSION_IDENTIFIER="${gene1}_${gene2}_${confidence}"
    SELECTED_FUSION_FILE="$OUTPUT_DIR/selected_fusion_${FUSION_IDENTIFIER}.tsv"
    VISUALIZATION_FILE="$OUTPUT_DIR/${REF_TYPE}_${FUSION_IDENTIFIER}_visualized.pdf"

    # Create R script to extract the specific fusion
    cat > "$OUTPUT_DIR/extract_fusion.R" << 'EOF'
args <- commandArgs(trailingOnly=TRUE)
fusion_file <- args[1]
gene1 <- args[2]
gene2 <- args[3]
direction1 <- args[4]
direction2 <- args[5]
confidence <- args[6]
output_file <- args[7]

# Read fusion file and clean column names
fusions <- read.delim(fusion_file, stringsAsFactors=FALSE)
names(fusions) <- gsub("X\\.", "", names(fusions))
names(fusions) <- gsub("\\.+", "", names(fusions))
names(fusions) <- gsub("genefusion", "(gene/fusion)", names(fusions))

# Find matching fusion with cleaned column names
matches <- fusions[
    toupper(fusions$gene1) == toupper(gene1) &
    toupper(fusions$gene2) == toupper(gene2) &
    fusions$`strand1(gene/fusion)` == direction1 &
    fusions$`strand2(gene/fusion)` == direction2 &
    tolower(fusions$confidence) == tolower(confidence),
]

if(nrow(matches) == 0) {
    cat("No matching fusion found with these exact criteria\n")
    quit(status=1)
}

# Get the header line from the original file
header_line <- readLines(fusion_file, n=1)

# Write the header and the matched row to the output file
writeLines(header_line, output_file)
write.table(matches[1,], 
            file=output_file, 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE,
            append=TRUE)

# Print confirmation of extraction
cat(sprintf("Successfully extracted fusion between %s and %s\n", gene1, gene2))
cat(sprintf("Breakpoints: %s - %s\n", matches$breakpoint1[1], matches$breakpoint2[1]))
cat(sprintf("Type: %s\n", matches$type[1]))
cat(sprintf("Confidence: %s\n", matches$confidence[1]))
EOF

    # Extract and visualize the specific fusion
    if Rscript "$OUTPUT_DIR/extract_fusion.R" "$FUSION_FILE" "$gene1" "$gene2" "$direction1" "$direction2" "$confidence" "$SELECTED_FUSION_FILE"; then
        echo "Generating visualization for selected fusion..."
        Rscript $ARRIBA_DIR/draw_fusions.R \
            --fusions="$SELECTED_FUSION_FILE" \
            --alignments="$BAM_FILE" \
            --output="$VISUALIZATION_FILE" \
            --annotation="$GTF_FILE" \
            --cytobands="$ARRIBA_DIR/database/cytobands_hg38_GRCh38_v2.4.0.tsv" \
            --proteinDomains="$ARRIBA_DIR/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
        echo "Visualization complete! Results saved in: $VISUALIZATION_FILE"
    else
        echo "No fusion found matching your criteria. Please check the fusion list and try again."
    fi

    # Ask if user wants to visualize another fusion
    read -p "Would you like to visualize another fusion? (y/n): " continue_answer
    if [[ ${continue_answer,,} != "y" ]]; then
        break
    fi
done

    # Update version info
    VERSION_INFO+="Fusion Visualization Only (${REF_TYPE}) - ${PATIENT_NAME}"
    echo "$VERSION_INFO" >> "$PATIENT_BASE_DIR/version_history.txt"
    exit 0
fi

# Create logs directory and input tracking file
mkdir -p "$OUTPUT_DIR/logs"
LOG_FILE="$OUTPUT_DIR/logs/analysis_$(date +%Y%m%d_%H%M%S).log"

# Check for FASTQ files
PATIENT_DIR="${BASE_DIR}/${PATIENT_NAME}"
R1_FILE=$(find "$PATIENT_DIR" -type f -name "*_*1.fq*.gz" | head -n 1)
R2_FILE=$(find "$PATIENT_DIR" -type f -name "*_*2.fq*.gz" | head -n 1)

if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ]; then
    echo "Error: Could not find paired FASTQ files in ${PATIENT_DIR}"
    exit 1
fi

# Create input tracking file
cat > "$OUTPUT_DIR/logs/input_info.txt" << EOF
=== Analysis Information ===
Date: $(date)
Patient Identifier: ${PATIENT_NAME}
Version: v${VERSION}

Input Files Used:
Forward Read (R1): $(basename "$R1_FILE")
Full Path: $R1_FILE

Reverse Read (R2): $(basename "$R2_FILE")
Full Path: $R2_FILE
EOF

echo "Found FASTQ files:"
echo "Forward reads: $(basename "$R1_FILE")"
echo "Reverse reads: $(basename "$R2_FILE")"

FASTQ_1="$R1_FILE"
FASTQ_2="$R2_FILE"

# Initialize version info
VERSION_INFO="v${VERSION} - $(date '+%Y-%m-%d %H:%M') - "

echo "Which analysis would you like to perform?"
echo "1) Fusion Analysis"
echo "2) Expression Analysis"
read -p "Enter your choice (1 or 2): " choice
validate_input "$choice" 1 2

# Only ask for reference choice if doing fusion analysis
if [ "$choice" = "1" ]; then
    echo "Which reference would you like to use?"
    echo "1) Ensembl (optimized for cancer fusions and internal tandem duplications like FLT3-ITD)"
    echo "2) Gencode (optimized for immune receptor gene fusions - Ig/TcR analysis)"
    read -p "Enter your choice (1 or 2): " ref_choice
    validate_input "$ref_choice" 1 2

   if [ "$ref_choice" = "1" ]; then
    GENOME_DIR=$ENSEMBL_DIR
    GENOME_FASTA=$ENSEMBL_FASTA
    GTF_FILE=$ENSEMBL_GTF
    REF_TYPE="ensembl"
    log_and_update_history "Selected Ensembl reference for cancer fusion and ITD analysis" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
    elif [ "$ref_choice" = "2" ]; then
    GENOME_DIR=$GENCODE_DIR
    GENOME_FASTA=$GENCODE_FASTA
    GTF_FILE=$GENCODE_GTF
    REF_TYPE="gencode"
    log_and_update_history "Selected Gencode reference for Ig/TcR analysis" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
    else
        echo "Invalid reference choice"
        exit 1
    fi
else
    # For expression analysis, set a default REF_TYPE as it's not used
    REF_TYPE="default"
fi

if [ "$choice" = "1" ]; then
    VERSION_INFO+="Fusion Analysis - ${REF_TYPE}"
    echo "Running Fusion Analysis for $PATIENT_NAME..."
  
  # Ask about QC/trimming
    echo -e "\nDo you want to perform quality control and trimming?"
    echo "Only choose yes if you have identified quality issues with FastQC"
    echo "1) No - proceed with original reads (recommended for high-quality samples)"
    echo "2) Yes - perform adapter trimming and quality filtering"
    read -p "Enter your choice (1 or 2): " qc_choice
    validate_input "$qc_choice" 1 2

    # Set up input files based on QC choice
    if [ "$qc_choice" = "2" ]; then
        echo "Running QC pipeline..."
        log_and_update_history "Starting QC pipeline" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
        run_qc_pipeline "$FASTQ_1" "$FASTQ_2" "$OUTPUT_DIR"
        log_and_update_history "Completed QC pipeline" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
        input_r1="${OUTPUT_DIR}/qc/cleaned_R1.fastq.gz"
        input_r2="${OUTPUT_DIR}/qc/cleaned_R2.fastq.gz"
    # Add QC information to version info
        VERSION_INFO+="-with_QC"
        
        # Log QC decision
        echo -e "\nQuality Control Applied:" >> "$OUTPUT_DIR/logs/input_info.txt"
        echo "Original R1: $(basename "$FASTQ_1")" >> "$OUTPUT_DIR/logs/input_info.txt"
        echo "Original R2: $(basename "$FASTQ_2")" >> "$OUTPUT_DIR/logs/input_info.txt"
        echo "QC-processed files used for analysis" >> "$OUTPUT_DIR/logs/input_info.txt"
    else
        echo "Proceeding with original reads..."
        input_r1="$FASTQ_1"
        input_r2="$FASTQ_2"
        
        # Log QC decision
        echo -e "\nQuality Control: Skipped" >> "$OUTPUT_DIR/logs/input_info.txt"
        echo "Using original FASTQ files" >> "$OUTPUT_DIR/logs/input_info.txt"
    fi


    log_and_update_history "Starting STAR alignment" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
  $STAR_DIR/STAR --runThreadN 8 \
    --genomeDir $GENOME_DIR \
    --readFilesIn "$input_r1" "$input_r2" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $OUTPUT_DIR/${REF_TYPE}_fusionaligned. \
    --readFilesCommand zcat \
    --outReadsUnmapped None \
    --outFilterMultimapNmax 10 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50 \
    --outBAMcompression 0 \
    --outFilterScoreMinOverLread 0.3 \
    --outFilterMatchNminOverLread 0.3 \
    --alignIntronMax 100000 \
    --alignMatesGapMax 100000 

# son 2yi sonradan ekledim, outfiltermultimapnmaxı da 50den 10a düşürdüm iyi çalışıyo dursun

    log_and_update_history "Starting Arriba fusion detection" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
     
    # Always run Arriba with blacklist first
$ARRIBA_DIR/arriba \
    -x $OUTPUT_DIR/${REF_TYPE}_fusionaligned.Aligned.sortedByCoord.out.bam \
    -o $OUTPUT_DIR/${REF_TYPE}_fusions_blacklist.tsv \
    -a $GENOME_FASTA \
    -g $GTF_FILE \
    -b $ARRIBA_DIR/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \
    -k $ARRIBA_DIR/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
    -t $ARRIBA_DIR/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
    -p $ARRIBA_DIR/database/protein_domains_hg38_GRCh38_v2.4.0.gff3

# Create interactive table for blacklist results
log_and_update_history "Creating interactive table for blacklisted fusions" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
Rscript create_fusion_table.R \
    "$OUTPUT_DIR/${REF_TYPE}_fusions_blacklist.tsv" \
    "$OUTPUT_DIR"
    
# Only if QC is chosen (contamination suspected), run additional no-blacklist analysis
if [ "$qc_choice" = "2" ]; then
    echo "Running additional no-blacklist analysis with specialized parameters..."
    # Prepare custom fusion database
    echo "Preparing custom fusion database..."
    zcat "$ARRIBA_DIR/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz" > "$OUTPUT_DIR/combined_fusions.tsv"
    cat "/home/viagen/rnaseq/arriba_v2.4.0/database/known_fusions_Gokce.tsv" >> "$OUTPUT_DIR/combined_fusions.tsv"

    $ARRIBA_DIR/arriba \
        -x $OUTPUT_DIR/${REF_TYPE}_fusionaligned.Aligned.sortedByCoord.out.bam \
        -o $OUTPUT_DIR/${REF_TYPE}_fusions_noblacklist.tsv \
        -a "$GENOME_FASTA" \
        -g "$GTF_FILE" \
        -k "$OUTPUT_DIR/combined_fusions.tsv" \
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

    # Create interactive table for no-blacklist results
    log_and_update_history "Creating interactive table for non-blacklisted fusions" "$LOG_FILE" "$PATIENT_BASE_DIR/version_history.txt"
    Rscript create_fusion_table.R \
        "$OUTPUT_DIR/${REF_TYPE}_fusions_noblacklist.tsv" \
        "$OUTPUT_DIR"

fi

# Always generate visualization for blacklisted results
echo "Generating visualization for blacklisted results..."
samtools index "$OUTPUT_DIR/${REF_TYPE}_fusionaligned.Aligned.sortedByCoord.out.bam"

Rscript $ARRIBA_DIR/draw_fusions.R \
    --fusions="$OUTPUT_DIR/${REF_TYPE}_fusions_blacklist.tsv" \
    --alignments="$OUTPUT_DIR/${REF_TYPE}_fusionaligned.Aligned.sortedByCoord.out.bam" \
    --output="$OUTPUT_DIR/${REF_TYPE}_visualizedfusions_blacklist.pdf" \
    --annotation="$GTF_FILE" \
    --cytobands="$ARRIBA_DIR/database/cytobands_hg38_GRCh38_v2.4.0.tsv" \
    --proteinDomains="$ARRIBA_DIR/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"

elif [ "$choice" = "2" ]; then
    echo "Running Expression Analysis for $PATIENT_NAME..."

    # Check for ReadsPerGene file
    READSPERGENE_FILE=$(find "${BASE_DIR}/${PATIENT_NAME}" -type f -iname "*ReadsPerGene.out.tab" | head -1)
    if [ ! -f "$READSPERGENE_FILE" ]; then
        echo "Error: Could not find ReadsPerGene.out.tab file in ${BASE_DIR}/${PATIENT_NAME}"
        exit 1
    fi
    READSPERGENE_FILENAME=$(basename "$READSPERGENE_FILE")
    echo -e "\nExpression Data File Used:" >> "$OUTPUT_DIR/logs/input_info.txt"
    echo "ReadsPerGene File: $READSPERGENE_FILE" >> "$OUTPUT_DIR/logs/input_info.txt"

    echo "What type of expression analysis would you like to perform?"
    echo "1) Compare with healthy tissue (using GTEx reference)"
    echo "2) Analyze tumor characteristics (using TCGA reference)"
    read -p "Enter your choice (1 or 2): " exp_choice
    validate_input "$exp_choice" 1 2

    if [ "$exp_choice" = "1" ]; then
         # Run GTEx analysis
    source gtex_analysis.sh "$PATIENT_NAME" "$READSPERGENE_FILE" "$OUTPUT_DIR" "$REF_TYPE" "$VERSION_INFO" "$READSPERGENE_FILENAME"
    if [ $? -ne 0 ]; then
        echo "GTEx analysis failed. Check the error messages above."
        exit 1
    fi
    elif [ "$exp_choice" = "2" ]; then
        # Run TCGA analysis - MODIFIED SECTION
        source tcga_analysis_fallback2025.sh "$PATIENT_NAME" "$READSPERGENE_FILE" "$OUTPUT_DIR" "$REF_TYPE" "$VERSION_INFO" "$READSPERGENE_FILENAME"
        # Note: Changed from VERSION_INFO=$(...) to source because we need the environment variables
        
        if [ $? -ne 0 ]; then
            echo "TCGA analysis failed. Check the error messages above."
            exit 1
        fi
    fi

# Run Reactome analysis after getting reference data
    if [ -f "$OUTPUT_DIR/normal_reference_counts.tab" ]; then
        echo "Running Reactome Analysis..."
        Rscript reactome_analysis_standalone_Haiku.R \
            "$OUTPUT_DIR/normal_reference_counts.tab" \
            "$READSPERGENE_FILE" \
            "$PATIENT_NAME"

        # Run KEGG pathway analysis if significant genes were found
if [ -f "$OUTPUT_DIR/expression_results/significant_genes.tsv" ]; then
    echo "Found significant_genes.tsv at: $OUTPUT_DIR/expression_results/significant_genes.tsv"
    echo "Running KEGG Pathway Analysis..."
    Rscript genelist3.r "$OUTPUT_DIR/expression_results/significant_genes.tsv" "$OUTPUT_DIR"
else
    echo "Error: Could not find significant_genes.tsv at: $OUTPUT_DIR/expression_results/significant_genes.tsv"
fi
    fi

# Add cleanup at the end of the script
rm -f "$OUTPUT_DIR/tcga_selections.rds"
rm -f "$OUTPUT_DIR/filtered_samples.rds"

    # Add version info to history file here
    echo "$VERSION_INFO" >> "$PATIENT_BASE_DIR/version_history.txt"

fi

