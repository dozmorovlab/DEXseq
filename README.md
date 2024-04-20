# Differential Exon Usage Analysis with DEXSeq

[DEXSeq Package](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)


## Analysis

[`00_Settings.sh`](./00_Settings.sh): Various settings for pipeline
- Dependencies
    - create a conda environment with:
        ```
        conda create --name DEU
        conda activate DEU
        conda install sra-tools samtools star r-base
        ```
- Parameters
    - `accesions`: a list of accessions to pull
    - `RNA_SEQ_DIR`: This is where a directory containing RNA-seq FASTQ files is located
    - `OUT_DIR`: This is where the output directory should be located
    - `MAX_CPUS`: Number of threads to use
    - `MAX_MEM_GB`: Maximum memory in GB to use
    - `HUMAN_GENOME_REF`: path to `.fna` FASTA reference file
    - `HUMAN_GENOME_GTF`: path to the GTF annotation file
    - `HUMAN_RNA_REF`: unused
    - `sjdbOverhang`: Set to the RNA read-length - 1 [more info](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

[`01_Download_Accessions.sh`](./01_Download_Accessions.sh): A script to download a given set of accessions with [sra-tools](https://github.com/ncbi/sra-tools)
- Output
    - FASTQs from given accesssions

[`02_Alignment.sh`](./02_Alignment.sh): A script to align the resulting FASTQs with STAR
- Output
    - BAMs for each sample; indexed and sorted

[`DEU_1.R`](./DEU.R): A script to perform Differential Exon Usage analysis with [DEXSeq](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)
- Input
    - BAMs
    - GTF for reference
- Output
    - `DEXSeqReport`: Standard `DEXSeq` report of all genes and their Differential Exon Usage (DEU) [more info](https://bioconductor.org/packages/release/bioc/manuals/DEXSeq/man/DEXSeq.pdf)
    - `dxr.rds`: a R dataset object containing the `DEXSeqResults` object generated during analysis

[`DEU_2.R`](./DEU.R): A script to analyze the resulting `DEXSeqResults` `.rds` file
- Input
    - `DEXSeqResults`: object generated during previous analysis
- Output
    - `DEU.xlsx`: Contains two sheets:
        - `Combined Genes`: For each gene, lists the following:
            - `gene`: gene name / alias
            - `gene_desc`: short description of the gene 
            - `total_exons`: total number of exons in the gene
            - `pvalue_sig_exons < 0.05`: number of exons with p-values that are significant in that gene given the threshold
            - `pvalue_prop_sig`: a proportion; $\frac{total\_exons}{pvalue\_sig\_exons < 0.05}$
            - `pvalue_na`: number of exons without a `pvalue`
            - `comb_pvalue`: the combined exon `pvalue` after using [Fisher's method](https://doi.org/10.2307%2F2681650)
            - `comb_pvalue_df`: degrees of freedom when combining each exon `pvalue`
            - `padj_sig_exons < 0.05`: number of exons with adjusted p-values that are significant in that gene given the threshold
            - `padj_prop_sig`: a proportion; $\frac{total\_exons}{padj\_sig\_exons < 0.05}$
            - `padj_na`: number of exons without a `padj`
            - `comb_padj`: the combined exon `padj` after using [Fisher's method](https://doi.org/10.2307%2F2681650)
            - `comb_padj_df`: degrees of freedom when combining each exon `padj`
        - `All Genes`: For each exon in each gene, lists the followning:
            - `groupID`: gene name / alias
            - `gene_desc`: short description of the gene 
            - `featureID`: feature ID
            - `exonBaseMean`: exon base mean
            - `dispersion`: dispersion
            - `stat`: statistic
            - `pvalue`: p-value
            - `padj`: adjusted p-value
            - `day_0`: condition #1
            - `day_100`: condition #2
            - `log2fold_day_100_day_0`: the $log_2$ fold change of the conditions; in this case $\frac{day\_100}{day\_0}$