library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(DEXSeq)

txdb <- makeTxDbFromGFF(
  "~/Ancestry/References/Human/ncbi_dataset/data/GCF_000001405.40/genomic.gtf"
)

flattenedAnnotation <- exonicParts(txdb, linked.to.single.gene.only = TRUE)
names(flattenedAnnotation) <- sprintf(
  "%s:E%0.3d",
  flattenedAnnotation$gene_id,
  flattenedAnnotation$exonic_part
)

# Count reads for each BAM
bamFiles_v <- c(
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435654.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435653.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435642.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435629.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435628.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435627.Aligned.sortedByCoord.out.bam",
  "~/../juicer/DEU/retinal_organoids/Data/star/BAMs/SRR15435652.Aligned.sortedByCoord.out.bam"
)
bamFiles <- BamFileList(bamFiles_v)
# seqlevelsStyle(flattenedAnnotation) <- "UCSC"
se <- summarizeOverlaps(
  flattenedAnnotation,
  bamFiles,
  singleEnd = FALSE,
  fragments = TRUE,
  ignore.strand = TRUE
)

# Building the DEXSeqDataSet
conditions <- c(
  "day_0", "day_0", "day_0",
  "day_100", "day_100", "day_100", "day_100"
)
colData(se)$condition <- factor(conditions)
colData(se)$libType <- factor(c(
  "paired-end", "paired-end", "paired-end",
  "paired-end", "paired-end", "paired-end", "paired-end"
))
dxd <- DEXSeqDataSetFromSE(se, design = ~ sample + exon + condition:exon)

# Run standard analysis
dxr <- DEXSeq(
  dxd,
  BPPARAM = MulticoreParam(workers = 64)
)
saveRDS(dxr, file = "~/../juicer/DEU/dxr.rds")

dxr <- readRDS("./dxr.RDS")

# Significant results
table(dxr$padj < 0.1)

# Order the subsetted data frame by padj in ascending order
dxr_by_padj <- dxr[order(dxr$padj), ]

length(unique(dxr_by_padj$groupID))

# Subset by specific rank
plotDEXSeq(
  dxr,
  unique(dxr_by_padj$groupID)[[5]],
  legend = TRUE,
  #cex.axis = 1,
  #cex = 1,
  #lwd = 1,
  #displayTranscripts = TRUE,
  names = TRUE,
  norCounts = TRUE,
  FDR = 0.1
)

DEXSeqHTML(
  dxr,
  FDR = 0.1,
  color = c("#FF000080", "#0000FF80"),
  BPPARAM = MulticoreParam(workers = 64)
)
