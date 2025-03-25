##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Support        : Code is based on following instructions: https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
## Date started 	: 2022-07-27
## Last modified  : 2025-03-24
## Aim 		        : Codes used for plotting Venn diagram
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.5D and Suppl.Fig.3A
##---------------------------------------------------------------------------

## Install ggvenn
#devtools::install_github("yanlinlin82/ggvenn")

## Load installed libraries
library(devtools)
library(ggplot2)
library(ggvenn)

## Empty workspace
rm(list=ls())

## Load data sets, For generation, please follow PPDPF_1_DeSeq2.R code
KO_day0 <- read.delim("DEGs_comparison_day0_KOvsWT.txt", header=T, sep="\t")
KO_day4 <- read.delim("DEGs_comparison_day4_KOvsWT.txt", header=T, sep="\t")
KO_day10 <- read.delim("DEGs_comparison_day10_KOvsWT.txt", header=T, sep="\t")
KO_day14 <- read.delim("DEGs_comparison_day14_KOvsWT.txt", header=T, sep="\t")

## Separate in upregulated and downregulated gene lists
KO_day0$revFC <- KO_day0$log2FoldChange * -1
KO_day0_up <-  KO_day0[which(KO_day0$padj<0.01 & KO_day0$log2FoldChange>1),] 
KO_day0_down <-  KO_day0[which(KO_day0$padj<0.01 & KO_day0$revFC>1),]

KO_day4$revFC <- KO_day4$log2FoldChange * -1
KO_day4_up <-  KO_day4[which(KO_day4$padj<0.01 & KO_day4$log2FoldChange>1),] 
KO_day4_down <-  KO_day4[which(KO_day4$padj<0.01 & KO_day4$revFC>1),]

KO_day10$revFC <- KO_day10$log2FoldChange * -1
KO_day10_up <-  KO_day10[which(KO_day10$padj<0.01 & KO_day10$log2FoldChange>1),] 
KO_day10_down <-  KO_day10[which(KO_day10$padj<0.01 & KO_day10$revFC>1),]

KO_day14$revFC <- KO_day14$log2FoldChange * -1
KO_day14_up <-  KO_day14[which(KO_day14$padj<0.01 & KO_day14$log2FoldChange>1),] 
KO_day14_down <-  KO_day14[which(KO_day14$padj<0.01 & KO_day14$revFC>1),]

## Create master list with all data to be plotted in Venn diagram
List_up <- list("0up" = KO_day0_up$X, "4up" =KO_day4_up$X, "10up" =KO_day10_up$X, "14up" =KO_day14_up$X)
List_down <- list("0down" = KO_day0_down$X, "4down" =KO_day4_down$X, "10down" =KO_day10_down$X, "14down" =KO_day14_down$X)

## Plot venn diagram
ggvenn(data = List_up,
  fill_color = c("#868686FF", "#41F2AF", "#7FE046", "#2A6304"),
  show_percentage = TRUE,digits = 0, fill_alpha = 0.75,stroke_size = 0.5, set_name_size = 4
)

ggvenn(data = List_down,
       fill_color = c("#868686FF", "#398DEA", "#1DC6ED", "#3131F7"),
       show_percentage = TRUE,digits = 0, fill_alpha = 0.75,stroke_size = 0.5, set_name_size = 4
)

###----------------------- part 2
## The following code was used for Suppl.Fig.3A to check intersection of DEGs with potential guide off-targets
## Check overlap of different nicking guides for a potential double strand break

guide1a <- read.delim("Guide1a_Potential_Off-target_cleaned.txt", header=F, sep="\t")
guide1b <- read.delim("Guide1b_Potential_Off-target_cleaned.txt", header=F, sep="\t")
guide2a <- read.delim("Guide2a_Potential_Off-target_cleaned.txt", header=F, sep="\t")
guide2b <- read.delim("Guide2b_Potential_Off-target_cleaned.txt", header=F, sep="\t")

List_guides <- list("1a" = guide1a$V1, "1b" =guide1b$V1, "2a" =guide2a$V1, "2b" =guide2b$V1)

## plot venn diagram
ggvenn(data = List_guides,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       show_percentage = TRUE,digits = 0, fill_alpha = 0.75,stroke_size = 0.5, set_name_size = 4
)

## Overlap and export of gene lists
overlap1 <- intersect(guide1a$V1, guide1b$V1)
# write.table(overlap1, "guide1ab_overlap.txt", quote=F, col.names=F, row.names=F)
overlap2 <- intersect(guide2a$V1, guide2b$V1)
# write.table(overlap2, "guide2ab_overlap.txt", quote=F, col.names=F, row.names=F)
All_guides <- union(
  union(guide1a$V1, guide1b$V1),
  union(guide2a$V1, guide2b$V1)
 )

## Additional comparisons that were considered but not depicted in the paper
# overlap1a2a <- intersect(guide1a$V1, guide2a$V1)
# write.table(overlap1a2a, "guide1a2a_overlap.txt", quote=F, col.names=F, row.names=F)
# overlap1a2b <- intersect(guide1a$V1, guide2b$V1)
# write.table(overlap1a2b, "guide1a2b_overlap.txt", quote=F, col.names=F, row.names=F)
# overlap1b2a <- intersect(guide1b$V1, guide2a$V1)
# write.table(overlap1b2a, "guide1b2a_overlap.txt", quote=F, col.names=F, row.names=F)
# overlap1b2b <- intersect(guide1b$V1, guide2b$V1)
# write.table(overlap1b2b, "guide1b2b_overlap.txt", quote=F, col.names=F, row.names=F)

###----------------------- part 3
## The following code was used for Suppl.Fig.3B to check intersection of genes that are DEGs at several time points with potential in silico predicted guide off-targets
## Check if any possible off-target genes overlap with up or downregulated genes
## Systematic comparison for at least one comp up, at least two, three or all...

## Overlap of all upregulated
## In all 4 stages up
day0_4_10_14_up <- intersect(
  intersect(KO_day0_up$X,KO_day4_up$X), 
  intersect(KO_day10_up$X,KO_day14_up$X))
## In at least 3 stages  up
up1 <- intersect(
  KO_day0_up$X, 
  intersect(KO_day4_up$X,KO_day10_up$X))
up2 <- intersect(
  KO_day0_up$X, 
  intersect(KO_day4_up$X,KO_day14_up$X))
up3 <- intersect(
  KO_day0_up$X, 
  intersect(KO_day10_up$X,KO_day14_up$X))
up4 <- intersect(
  KO_day4_up$X, 
  intersect(KO_day10_up$X,KO_day14_up$X))  
threedays_up <- union(
  union(up1, up2),
  union(up3, up4))
## In at least 2 stages up
Two_up1 <- intersect(KO_day0_up$X,KO_day4_up$X)
Two_up2 <- intersect(KO_day0_up$X,KO_day10_up$X)
Two_up3 <- intersect(KO_day0_up$X,KO_day14_up$X)
Two_up4 <- intersect(KO_day4_up$X,KO_day10_up$X)
Two_up5 <- intersect(KO_day4_up$X,KO_day14_up$X)
Two_up6 <- intersect(KO_day10_up$X,KO_day14_up$X)
twodays_up <- union(
  union(union(Two_up1, Two_up2),union(Two_up3, Two_up4)),
  union(Two_up5, Two_up6))
## In at least one stage up
up_0_4 <- union(KO_day0_up$X, KO_day4_up$X)
up_10_14 <- union(KO_day10_up$X, KO_day14_up$X)
All_up <- union(up_0_4, up_10_14)

## Master list for upregulated comp                   
List_guides_DEG_up <- list("guides" = All_guides, "1comp" =All_up, "2comp" = twodays_up, "3comp"= threedays_up)

## Plot Venn, overlap to guide depicted in Suppl.Fig.3B left
ggvenn(data = List_guides_DEG_up,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
       show_percentage = TRUE,digits = 0, fill_alpha = 0.75,stroke_size = 0.5, set_name_size = 4
)

## Overlap of all downregulated
## In all 4 stages down
day0_4_10_14_down <- intersect(
  intersect(KO_day0_down$X,KO_day4_down$X), 
  intersect(KO_day10_down$X,KO_day14_down$X))
## In at least 3 stages  down
down1 <- intersect(
  KO_day0_down$X, 
  intersect(KO_day4_down$X,KO_day10_down$X))
down2 <- intersect(
  KO_day0_down$X, 
  intersect(KO_day4_down$X,KO_day14_down$X))
down3 <- intersect(
  KO_day0_down$X, 
  intersect(KO_day10_down$X,KO_day14_down$X))
down4 <- intersect(
  KO_day4_down$X, 
  intersect(KO_day10_down$X,KO_day14_down$X))  
threedays_down <- union(
  union(down1, down2),
  union(down3, down4))
## In at least 2 stages down
Two_down1 <- intersect(KO_day0_down$X,KO_day4_down$X)
Two_down2 <- intersect(KO_day0_down$X,KO_day10_down$X)
Two_down3 <- intersect(KO_day0_down$X,KO_day14_down$X)
Two_down4 <- intersect(KO_day4_down$X,KO_day10_down$X)
Two_down5 <- intersect(KO_day4_down$X,KO_day14_down$X)
Two_down6 <- intersect(KO_day10_down$X,KO_day14_down$X)
twodays_down <- union(
  union(union(Two_down1, Two_down2),union(Two_down3, Two_down4)),
  union(Two_down5, Two_down6))
## In at least one stage down
down_0_4 <- union(KO_day0_down$X, KO_day4_down$X)
down_10_14 <- union(KO_day10_down$X, KO_day14_down$X)
All_down <- union(down_0_4, down_10_14)

## Master list for downregulated comp                   
List_guides_DEG_down <- list("guides" = All_guides, "1comp" =All_down, "2comp" = twodays_down, "3comp"= threedays_down)

## Plot Venn, overlap to guide depicted in Sdownpl.Fig.3B left
ggvenn(data = List_guides_DEG_down,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
       show_percentage = TRUE,digits = 0, fill_alpha = 0.75,stroke_size = 0.5, set_name_size = 4
)
