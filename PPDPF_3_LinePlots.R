##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Support        : Coding Template from Kevin Menden (ORCID: 0000-0002-5228-985X)
## Date started 	: 2019-10-06
## Last modified  : 2025-03-24
## Aim 		        : Codes used for gene expression line plots from preprocessed and normalized RNA-seq bulk data
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.6E-F and Suppl.Fig.3E,4)
##--------------------------------------------------------------------------- 

## Load installed libraries
library(ggplot2)
library(dplyr)

## Empty workspace
rm(list=ls())

## Load normalized expression data
mat <- read.table("GSE282625_ReadCounts_deseq-Rlognormalized.txt", sep="\t", header = T, row.names = 1)
legend<- read.delim("legend.txt", sep="\t", header = T)

## 1. Loop for plotting and saving line graphs for all genes
iter <- length(rownames(mat))
for (i in 1:iter) {
  print(i)
  gene_of_interest <- rownames(mat)[i]
  ## Subset for gene
  gene_exp <- mat[gene_of_interest,]
  ## Divide WT and KO
  group <- rep("Control", ncol(gene_exp))
  group[grepl("KO", colnames(gene_exp))] <- "KO"
  group_clone <- legend$clone
  
  ## Vector for the days
  day <- colnames(gene_exp)
  day[grepl("day4", colnames(gene_exp))] <- 4
  day[grepl("day0", colnames(gene_exp))] <- 0
  day[grepl("day10", colnames(gene_exp))] <- 10
  day[grepl("day14", colnames(gene_exp))] <- 14
  day <- as.numeric(day)
  
  df <- data.frame(exp = as.numeric(gene_exp), group = group_clone, day = day)
  df <- df[order(df$day),]
  
  ## Summarize mean
  agg_df <- aggregate.data.frame(df, by=list(df$day, df$group), FUN = mean)
  agg_df <- agg_df[order(agg_df$day),]
  agg_df$group <- agg_df$Group.2
  
  ## Add standard Error
  ste_mean <- function(x) { sd(x)/sqrt(length(x))}
  agg_df_stE <- aggregate.data.frame(df, by=list(df$day, df$group), FUN = ste_mean)
  agg_df_stE <- agg_df_stE[order(agg_df_stE$day),]
  agg_df_stE$group <- agg_df_stE$Group.2
  agg_df$stE <- agg_df_stE$exp
  
  p <- ggplot(data = agg_df, aes(x=day, y=exp, group=group, color=group)) +
    geom_line(linetype="twodash") +
    theme_classic() +
    scale_color_manual(values=c('#7FE046','#3D8E09', '#1DC6ED','#398DEA')) + 
    geom_errorbar(aes(ymin=exp-stE, ymax=exp+stE), width=.2) + 
    ggtitle(gene_of_interest) +
    scale_x_discrete(name ="Day of differentiation", limits=c(0, 4, 10, 14)) +
    ylab("Normalized expression (R-log)")
  
  ggsave(filename = paste("plots/",gene_of_interest, "_plot.pdf", sep=""), plot = p)
}

## 2. Alternative to Loop: Manual plotting for selected genes
gene_of_interest <- "TRPC7"
## Subset for gene
gene_exp <- mat[gene_of_interest,]
## Divide WT and KO
group <- rep("Control", ncol(gene_exp))
group[grepl("KO", colnames(gene_exp))] <- "KO"
group_clone <- legend$clone

## Vector for the days
day <- colnames(gene_exp)
day[grepl("day4", colnames(gene_exp))] <- 4
day[grepl("day0", colnames(gene_exp))] <- 0
day[grepl("day10", colnames(gene_exp))] <- 10
day[grepl("day14", colnames(gene_exp))] <- 14
day <- as.numeric(day)

df <- data.frame(exp = as.numeric(gene_exp), group = group_clone, day = day)
df <- df[order(df$day),]

## Summarize mean
agg_df <- aggregate.data.frame(df, by=list(df$day, df$group), FUN = mean)
agg_df <- agg_df[order(agg_df$day),]
agg_df$group <- agg_df$Group.2

## Add standard Error
ste_mean <- function(x) { sd(x)/sqrt(length(x))}
agg_df_stE <- aggregate.data.frame(df, by=list(df$day, df$group), FUN = ste_mean)
agg_df_stE <- agg_df_stE[order(agg_df_stE$day),]
agg_df_stE$group <- agg_df_stE$Group.2
agg_df$stE <- agg_df_stE$exp

p <- ggplot(data = agg_df, aes(x=day, y=exp, group=group, color=group)) +
  geom_line(linetype="twodash") +
  theme_classic() +
  scale_color_manual(values=c('#7FE046','#3D8E09', '#1DC6ED','#398DEA')) + 
  geom_errorbar(aes(ymin=exp-stE, ymax=exp+stE), width=.2) + 
  ggtitle(gene_of_interest) +
  scale_x_discrete(name ="Day of differentiation", limits=c(0, 4, 10, 14)) +
  ylab("Normalized expression (R-log)")
p



