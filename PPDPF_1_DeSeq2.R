##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Support        : Coding Templates from Kevin Menden (ORCID: 0000-0002-5228-985X)
## Date started 	: 2019-10-06
## Last modified  : 2025-03-24
## Aim 		        : Codes used for analyzing pre-processed bulk RNA-seq data
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Supplementary Table 3 and foundation for other codes (Fig.5-7 and Suppl.Fig.3-4)
##--------------------------------------------------------------------------- 

## Load installed libraries
library(DESeq2)
library(biomaRt)

## Empty workspace
rm(list = ls())

## Load data
counts <- read.table("GSE282625_ReadCounts_NotNormalized.txt", sep="\t", header = T, row.names = 1)

## Set Cutoff values
row_sum_cutoff = 30
pvalue_cutoff = 0.05
lfc_cutoff = 1

## Convert Ensembl IDs to HGNC gene names
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- rownames(counts)
bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values = genes, mart = ensembl)
mdata <- merge(counts, bm, by.x="row.names", by.y="ensembl_gene_id")
mdata <- mdata[!duplicated(mdata$hgnc_symbol),]
rownames(mdata) <- mdata$hgnc_symbol
new_counts <- mdata[,c(-1, -ncol(mdata))]
counts <- new_counts

## Define group names as input für DeSeq2
cols <- colnames(counts)
groups <- as.character(sapply(cols, FUN = function(x){strsplit(x, split="HUES8")[[1]][[2]]}))
groups <- gsub("_EXDPF_", "", groups)
groups <- make.names(groups)
groups <- gsub("X_", "", groups)
df <- data.frame(group = groups)
rownames(df) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(counts,
                              colData = df,
                              design = ~ group)
keep <- rowSums(counts(dds)) >= row_sum_cutoff
dds <- dds[keep,]
dds <- DESeq(dds)

## Do WT KO Comparisons for each stage
## Comparison day0
res_day0 <- results(dds, contrast = c("group", "KO_day0", "Control_day0"))
res_day0 <- na.omit(res_day0)
res_day0.sig <- res_day0[res_day0$padj <= 0.05,]
res_day0.sig <- res_day0.sig[abs(res_day0.sig$log2FoldChange) >= 1,]

## Comparison day4
res_day4 <- results(dds, contrast = c("group", "KO_day4", "Control_day4"))
res_day4 <- na.omit(res_day4)
res_day4.sig <- res_day4[res_day4$padj <= 0.05,]
res_day4.sig <- res_day4.sig[abs(res_day4.sig$log2FoldChange) >= 1,]

## Comparison day day 10
res_day10 <- results(dds, contrast = c("group", "KO_day10", "Control_day10"))
res_day10 <- na.omit(res_day10)
res_day10.sig <- res_day10[res_day10$padj <= 0.05,]
res_day10.sig <- res_day10.sig[abs(res_day10.sig$log2FoldChange) >= 1,]

## Comparison day day 14
res_day14 <- results(dds, contrast = c("group", "KO_day14", "Control_day14"))
res_day14 <- na.omit(res_day14)
res_day14.sig <- res_day14[res_day14$padj <= 0.05,]
res_day14.sig <- res_day14.sig[abs(res_day14.sig$log2FoldChange) >= 1,]

## Comparison of WT over time
## Comparison WT day4-day0
res_wt4to0 <- results(dds, contrast = c("group", "Control_day4", "Control_day0"))
res_wt4to0 <- na.omit(res_wt4to0)
res_wt4to0.sig <- res_wt4to0[res_wt4to0$padj <= 0.05,]
res_wt4to0.sig <- res_wt4to0.sig[abs(res_wt4to0.sig$log2FoldChange) >= 1,]

## Comparison WT day10-day0
res_wt10to0 <- results(dds, contrast = c("group", "Control_day10", "Control_day0"))
res_wt10to0 <- na.omit(res_wt10to0)
res_wt10to0.sig <- res_wt10to0[res_wt10to0$padj <= 0.05,]
res_wt10to0.sig <- res_wt10to0.sig[abs(res_wt10to0.sig$log2FoldChange) >= 1,]

## Comparison WT day14-day0
res_wt14to0 <- results(dds, contrast = c("group", "Control_day14", "Control_day0"))
res_wt14to0 <- na.omit(res_wt14to0)
res_wt14to0.sig <- res_wt14to0[res_wt14to0$padj <= 0.05,]
res_wt14to0.sig <- res_wt14to0.sig[abs(res_wt14to0.sig$log2FoldChange) >= 1,]

## Comparison WT day10-day4
res_wt10o4 <- results(dds, contrast = c("group", "Control_day10", "Control_day4"))
res_wt10o4 <- na.omit(res_wt10o4)
res_wt10o4.sig <- res_wt10o4[res_wt10o4$padj <= 0.05,]
res_wt10o4.sig <- res_wt10o4.sig[abs(res_wt10o4.sig$log2FoldChange) >= 1,]

## Comparison WT day14-day10
res_wt14o10 <- results(dds, contrast = c("group", "Control_day14", "Control_day10"))
res_wt14o10 <- na.omit(res_wt14o10)
res_wt14o10.sig <- res_wt14o10[res_wt14o10$padj <= 0.05,]
res_wt14o10.sig <- res_wt14o10.sig[abs(res_wt14o10.sig$log2FoldChange) >= 1,]

## Comparison WT day14-day4
res_wt14to4 <- results(dds, contrast = c("group", "Control_day14", "Control_day4"))
res_wt14to4 <- na.omit(res_wt14to4)
res_wt14to4.sig <- res_wt14to4[res_wt14to4$padj <= 0.05,]
res_wt14to4.sig <- res_wt14to4.sig[abs(res_wt14to4.sig$log2FoldChange) >= 1,]

## save expression list in specified subfolder
write.table(res_day0, "DEGs_comparison_day0_KOvsWT.txt", sep="\t", quote=F, col.names = NA)
write.table(res_day4, "DEGs_comparison_day4_KOvsWT.txt", sep="\t", quote=F, col.names = NA)
write.table(res_day10, "DEGs_comparison_day10_KOvsWT.txt", sep="\t", quote=F, col.names = NA)
write.table(res_day14, "DEGs_comparison_day14_KOvsWT.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt4to0, "DEGs_comparison_WT_day4_vs_day0.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt10to0, "DEGs_comparison_WT_day10_vs_day0.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt14to0, "DEGs_comparison_WT_day14_vs_day0.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt10o4, "DEGs_comparison_WT_day10_vs_day4.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt14o10, "DEGs_comparison_WT_day14_vs_day10.txt", sep="\t", quote=F, col.names = NA)
write.table(res_wt14to4, "DEGs_comparison_WT_day14_vs_day4.txt", sep="\t", quote=F, col.names = NA)

## Extract normalized counts
#norm.counts <- counts(dds, normalized=TRUE)
#write.table(norm.counts, "deseq_normalized_counts.txt", sep="\t", quote=F, col.names = NA)

## Calculate VST and Rlog normalized counts for downstream applications and visualizations
#mat_VST <- varianceStabilizingTransformation(dds, blind = FALSE)
#mat_VST <- assay(mat_VST)
#write.table(mat_VST, "deseq_VSTnormalized_counts.txt", sep="\t", quote=F, col.names = NA)

mat_Rlog <- rlog(dds, blind = FALSE)
mat_Rlog <- assay(mat_Rlog)
write.table(mat_Rlog, "GSE282625_ReadCounts_deseq-Rlognormalized.txt", sep="\t", quote=F, col.names = NA)

##Plotting QC plots
plotDispEsts(dds)
barplot(colSums(counts(dds)))
barplot(colSums(counts(dds, normalized=T)))
