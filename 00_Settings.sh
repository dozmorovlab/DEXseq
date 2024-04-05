# Various settings for later scripts

# Activate conda environment
# 
# Ex:
# conda create --name DEU
# conda activate DEU
# conda install sra-tools samtools star r-base
conda init
conda activate DEU

# RNA-seq directory
RNA_SEQ_DIR=~/../juicer/DEU/retinal_organoids

# Output data directory
OUT_DIR=~/../juicer/DEU/retinal_organoids

# Resources
# qrsh -pe smp 64 -l mem_free=200G
MAX_CPUS=64
MAX_MEM_GB=200

# Human reference
# https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,RNA_FASTA,GENOME_GTF
HUMAN_GENOME_REF=~/Ancestry/References/Human/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
HUMAN_GENOME_GTF=~/Ancestry/References/Human/ncbi_dataset/data/GCF_000001405.40/genomic.gtf
#HUMAN_RNA_REF=~/Ancestry/References/Human/ncbi_dataset/data/GCF_000001405.40/rna.fna

# Set to the RNA read-length - 1
sjdbOverhang=149

# Create Directory Structure
# ---------------------------------------------------

mkdir -p $OUT_DIR

# Data
data_dir=$OUT_DIR/Data
mkdir -p $data_dir

star_data="$data_dir/star"
mkdir -p $star_data
star_index_dir="$star_data/index"
mkdir -p $star_index_dir
star_bam_dir="$star_data/BAMs"
mkdir -p $star_bam_dir

# Reports
report_dir=$OUT_DIR/Reports
mkdir -p $report_dir

star_report="$report_dir/star"
mkdir -p $star_report

# Functions
# ---------------------------------------------------

process_fastq_pairs() {
    local IN_RNA_SEQ_DIR="$1"
    local R1=${2:-"_R1_"}
    local R2=${3:-"_R2_"}

    local FASTQ_PAIRS=()

    # Find all FASTQ.gz files in the specified directory
    FILES=$(find "$IN_RNA_SEQ_DIR" -type f -name "*.fastq.gz")

    # Iterate over each file
    for FILE in $FILES; do
        # Get the file name without extension
        FILENAME=$(basename "$FILE" .fastq.gz)
        DIRNAME=$(dirname "$FILE")

        # Check if the file name contains "R1"
        if [[ "$FILENAME" == *"$R1" ]]; then
            # Remove "R1" from the file name to get the corresponding "R2" file name
            R2_FILENAME=${FILENAME/"$R1"/"$R2"}
            
            # Check if the corresponding "R2" file exists
            R2_FILE="$DIRNAME/$R2_FILENAME.fastq.gz"
            if [ -f "$R2_FILE" ]; then
                # Add the pair to the list
                FASTQ_PAIRS+=("$FILE" "$R2_FILE")
            fi
        fi
    done

    # Return the array of FASTQ pairs
    echo "${FASTQ_PAIRS[@]}"
}

