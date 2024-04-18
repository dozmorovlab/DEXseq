# Brydon Wall

require(org.Hs.eg.db)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(DEXSeq)
library(AnnotationDbi)
library(metap)
library(openxlsx)

pvalue_cutoff <- 0.05
padj_cutoff <- 0.05

dxr <- readRDS("./dxr.RDS")

col_exclude <- c(
  "genomicData",
  "countData",
  "transcripts"
)

# Exclude cols
dxr <- dxr[, -which(names(dxr) %in% col_exclude)]

# Make data.frame
all_df <- as.data.frame(dxr)

# Get descriptions
desc_select <- select(
  org.Hs.eg.db,
  keys = unique(all_df$groupID),
  keytype = "SYMBOL",
  columns = "GENENAME"
)

# Remove duplicates, keep first
desc_select <- desc_select[!duplicated(desc_select$SYMBOL), ]

# Adjust colnames
colnames(desc_select) <- c("temp", "gene_desc")

# Account for underscores for descriptions
all_df["temp"] <- gsub("_.*", "", all_df$groupID)

# Merge and reorder columns
all_df <- merge(all_df, desc_select, by = "temp")
all_df["temp"] <- NULL
column_order <- c(
  "groupID", "gene_desc", "featureID", "exonBaseMean",
  "dispersion", "stat", "pvalue", "padj", "day_0",
  "day_100", "log2fold_day_100_day_0"
)
all_df <- all_df[, column_order]

# Sort
all_df <- all_df[order(all_df$featureID), ]
all_df <- all_df[order(all_df$groupID), ]

genes_df <- data.frame(
  "gene" = c(),
  "gene_desc" = c(),
  "total_exons" = c(),

  "pvalue_sig_exons" = c(),
  "pvalue_prop_sig" = c(),
  "pvalue_na" = c(),
  "comb_pvalue" = c(),
  "comb_pvalue_df" = c(),

  "padj_sig_exons" = c(),
  "padj_prop_sig" = c(),
  "padj_na" = c(),
  "comb_padj" = c(),
  "comb_padj_df" = c()
)

# Genes
genes <- unique(all_df$groupID)
gene <- "DACH1"
genes_df <- do.call(rbind, lapply(genes, function(gene) {

  gene_df <- all_df[all_df$groupID == gene, ]
  gene_desc <- gene_df$gene_desc[[1]]
  total_exons <- nrow(gene_df)

  # pvalue
  pvalue_sig_exons <- sum(
    gene_df$pvalue[!is.na(gene_df$pvalue)] < pvalue_cutoff
  )
  pvalue_prop_sig <- pvalue_sig_exons / total_exons
  pvalue_log <- sumlog(
    gene_df$pvalue[!is.na(gene_df$pvalue)]
  )
  comb_pvalue <- pvalue_log[["p"]]
  comb_pvalue_df <- pvalue_log[["df"]]
  pvalue_na <- sum(is.na(gene_df$pvalue))

  # padj
  padj_sig_exons <- sum(
    gene_df$padj[!is.na(gene_df$padj)] < padj_cutoff
  )
  padj_prop_sig <- padj_sig_exons / total_exons
  padj_log <- sumlog(
    gene_df$padj[!is.na(gene_df$padj)]
  )
  comb_padj <- padj_log[["p"]]
  comb_padj_df <- padj_log[["df"]]
  padj_na <- sum(is.na(gene_df$padj))

  data.frame(
    gene = gene,
    gene_desc = gene_desc,
    total_exons = total_exons,

    pvalue_sig_exons = pvalue_sig_exons,
    pvalue_prop_sig = pvalue_prop_sig,
    pvalue_na = pvalue_na,
    comb_pvalue = comb_pvalue,
    comb_pvalue_df = comb_pvalue_df,

    padj_sig_exons = padj_sig_exons,
    padj_prop_sig = padj_prop_sig,
    padj_na = padj_na,
    comb_padj = comb_padj,
    comb_padj_df = comb_padj_df
  )
}))

genes_df <- genes_df[order(genes_df$comb_padj), ]

# Rename cols to add cutoffs
names(genes_df)[names(genes_df) == "pvalue_sig_exons"] <- paste0(
  "pvalue_sig_exons < ", pvalue_cutoff
)
names(genes_df)[names(genes_df) == "padj_sig_exons"] <- paste0(
  "padj_sig_exons < ", padj_cutoff
)

# Create a new Excel workbook
wb <- createWorkbook(creator = "Brydon P. G. Wall")

# Add the first sheet to the workbook
addWorksheet(wb, "Combined Genes")
writeData(wb, "Combined Genes", genes_df)

# Add the second sheet to the workbook
addWorksheet(wb, "All Genes")
writeData(wb, "All Genes", all_df)

# Save the workbook to a file
saveWorkbook(wb, "DEU.xlsx", overwrite = TRUE)
