# Differential Exon Usage Analysis with DEXSeq

[DEXSeq Package](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)


## Analysis

[`00_Settings.sh`](./00_Settings.sh): Various settings for pipeline
- Input
    - todo
- Output
    - todo

[`01_Download_Accessions.sh`](./01_Download_Accessions.sh): A script to download a given set of accessions with [sra-tools](https://github.com/ncbi/sra-tools)
- Input
    - todo
- Output
    - todo

[`02_Alignment.sh`](./02_Alignment.sh): A script to align the resulting FASTQs with STAR
- Input
    - todo
- Output
    - todo

[`DEU_1.R`](./DEU.R): A script to perform Differential Exon Usage analysis with [DEXSeq](https://doi.org/doi:10.18129/B9.bioc.DEXSeq)
- Input
    - todo
- Output
    - todo

[`DEU_2.R`](./DEU.R): A script to analyze the resulting `DEXSeqResults` `.rds` file
- Input
    - todo
- Output
    - todo

## Results

`DEXSeqReport`: Standard `DEXSeq` report of all genes and their Differential Exon Usage (DEU)

`DEU.xlsx`: Contains two sheets:
- `Combined Genes`: For each gene, lists the following:
    - gene
    - gene_desc
    - total_exons
    - pvalue_sig_exons < 0.05
    - pvalue_prop_sig
    - pvalue_na
    - comb_pvalue
    - comb_pvalue_df
    - padj_sig_exons < 0.05
    - padj_prop_sig
    - padj_na
    - comb_padj
    - comb_padj_df
- `All Genes`: For each exon in each gene, lists the followning:
    - groupID
    - gene_desc
    - featureID
    - exonBaseMean
    - dispersion
    - stat
    - pvalue
    - padj
    - day_0
    - day_100
    - log2fold_day_100_day_0