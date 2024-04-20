qrsh -pe smp 32
cd ~/../juicer/retinal_organoids
conda activate DEU
source ./00_Stttings.sh

threads=$MAX_CPUS

for accession in $accessions; do
    fasterq-dump $accession \
        --split-3 \
        --threads $threads
    
    fastqs=(
        "./${accession}_1.fastq"
        "./${accession}_2.fastq"
        "./${accession}.fastq"
    )
    
    for fastq in "${fastqs[@]}"; do
        if [ -f "$fastq" ]; then
            bgzip "$fastq" \
                -@ $threads
        fi
    done
done
