qrsh -pe smp 32
cd ~/../juicer/retinal_organoids
conda activate DEU

threads=32

accessions=(
    SRR15435654
    SRR15435653
    SRR15435642
    SRR15435629
    SRR15435628
    SRR15435627
    SRR15435652
)

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
