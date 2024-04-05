# Differential Exon Usage

[10.18129/B9.bioc.DEXSeq](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)



[00_Settings.sh](./00_Settings.sh): Various settings for pipeline

[01_Download_Accessions.sh](./01_Download_Accessions.sh): A script to download a given set of accessions with [sra-tools](https://github.com/ncbi/sra-tools)

[02_Alignment.sh](./02_Alignment.sh): A script to align the resulting FASTQs with STAR

[DEU.R](./DEU.R): A script to perform Differential Exon Usage analysis with [DEXSeq](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)
