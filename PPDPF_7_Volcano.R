##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Support        : Code is based on following github tutorial: #https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
## Date started 	: 2024-02-14
## Last modified  : 2025-03-24
## Aim 		        : Codes used for plotting Venn diagram
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.5F
##---------------------------------------------------------------------------

## Load installed libraries
library(ggrepel)
library(dplyr)
library(ggbreak) #for breaking x/y axes

## Empty workspace
rm(list=ls())

## Load data sets, For generation, please follow PPDPF_1_DeSeq2.R code
de_0 <- read.delim("DEGs_comparison_day0_KOvsWT.txt", header=T, sep="\t")
de_4 <- read.delim("DEGs_comparison_day4_KOvsWT.txt", header=T, sep="\t")
de_0 <- read.delim("DEGs_comparison_day10_KOvsWT.txt", header=T, sep="\t")
de_14 <- read.delim("DEGs_comparison_day14_KOvsWT.txt", header=T, sep="\t")

## Add a column to specify up/down direction
de_0$diffexpressed <- "NO"
de_0$diffexpressed[de_0$log2FoldChange > 1 & de_0$padj < 0.01] <- "UP"
de_0$diffexpressed[de_0$log2FoldChange < -1 & de_0$padj < 0.01] <- "DOWN"

## Only select top 10% FC of DEGs for labeling
de_0_up <- de_0 %>%
  filter(diffexpressed == "UP") %>%
  arrange(desc(log2FoldChange))
Top10_up <- de_0_up$X[1:10]

de_0_down <- de_0 %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(log2FoldChange)
Top10_down <- de_0_down$X[1:10]

## Label points with genes names
## Create a new column "delabel" to de, that will contain the name of Top10 DEGs
de_0$delabel <- NA
de_0$delabel[which(de_0$X %in% Top10_down)] <- de_0$X[which(de_0$X %in% Top10_down)]
de_0$delabel[which(de_0$X %in% Top10_up)] <- de_0$X[which(de_0$X %in% Top10_up)]

## The basic scatter plot: x is "log2FoldChange", y is log10 p value
## If color for significant proteins is desired, col=diffexpressed
p <- ggplot(data=de_0, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_classic() + 
  #xlim(-10,10) +
  #ylim(NA,40) +
  #scale_x_cut(breaks=c(-10)) +
  #scale_y_break(breaks=c(40,41)) +
  geom_text_repel(position = position_dodge(width = NULL, preserve = "total"),  size=(3)) + # geom_text_repel() for some, geom_text() for all sig labels
  #geom_hline(yintercept=-log10(0.05), col="black") +
  scale_color_manual(values=c("#1DC6ED",  "#000000", "#7FE046")) 

## Adapt axis and title and appearance
mynamestheme <- theme(
  legend.title = element_text(family = "Arial", size = (12)),
  legend.text = element_text(family = "Arial", size = (10)),
  axis.title = element_text(family = "Arial", size = (12)),
  axis.text = element_text(family = "Arial", size = (10))
)

## Plot Volcano plot (colour for legend title)
print(p+ mynamestheme + labs(colour= "DEGs", y= "-log10 P adj", x= "log2FC"))

## Save the plot
# pdf("Volcano_d10_YCut.pdf")
# print(p+labs(colour= "DEGs", y= "-log10 P adje", x= "log2FC"))
# dev.off() 

