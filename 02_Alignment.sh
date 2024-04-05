#!/usr/bin/env bash
#$ -cwd
# shell for qsub to use:
#$ -S /bin/bash
#$ -N Alignment_03
#$ -e ./01_Alignment_errors.txt
#$ -o ./01_Alignment_output.txt
#$ -pe smp 64
#$ -l mem_free=200G

# Load settings
source ./00_Settings.sh

# Index once (check for existing index)
if [ ! "$(ls -A $star_index_dir)" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing human genome: $HUMAN_GENOME_REF"
    STAR \
		--runThreadN $MAX_CPUS \
		--runMode genomeGenerate \
		--genomeDir $star_index_dir \
		--genomeFastaFiles $HUMAN_GENOME_REF \
		--sjdbGTFfile $HUMAN_GENOME_GTF \
		--sjdbOverhang $sjdbOverhang
fi

# Align with STAR
# ENCODE standard options for long RNA-seq pipeline
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
FASTQ_PAIRS=($(process_fastq_pairs "$RNA_SEQ_DIR" "_1" "_2"))
MAX_MEM_BYTES="$(($MAX_MEM_GB * 1073741824))"

for ((i = 0; i < ${#FASTQ_PAIRS[@]}; i+=2)); do
	R1=${FASTQ_PAIRS[i]}
    R2=${FASTQ_PAIRS[i+1]}
	neutral_base=$(basename "$R1" _1.fastq.gz)

	echo "$(date '+%Y-%m-%d %H:%M:%S'): Aligning $R1 and $R2"

	STAR \
		--runThreadN $MAX_CPUS \
		--limitBAMsortRAM $MAX_MEM_BYTES \
		--genomeDir $star_index_dir \
		--readFilesIn $R1 $R2 \
		--outFileNamePrefix "$star_bam_dir/${neutral_base}." \
		--outSAMunmapped Within KeepPairs \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMcompression 10 \
		--outBAMsortingThreadN $MAX_CPUS \
		--readFilesCommand "bgzip -cd -@ $MAX_CPUS" \
		--outFilterType BySJout \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--twopassMode Basic
	
	cp "$star_bam_dir/${neutral_base}.Log.final.out" "$star_report/${neutral_base}.Log.final.out"
done

# Index BAMs
for file in "$star_bam_dir"/*.bam; do
	echo "$(date '+%Y-%m-%d %H:%M:%S'): Indexing $file"
	samtools index "$file" \
	--threads $MAX_CPUS
done
