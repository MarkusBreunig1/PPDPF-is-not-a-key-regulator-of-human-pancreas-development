##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Date started 	: 2022
## Last modified  : 2025-03-24
## Aim 		        : Codes used for Principal Component analysis (PCA) on normalized RNA-seq bulk data
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.5A,B
##---------------------------------------------------------------------------

## Load installed libraries
library(RColorBrewer)
library(dplyr)
library("ggbiplot")
library("factoextra")

## Empty workspace
rm(list=ls())

## Load normalized expression data and legend
counts <- read.delim("GSE282625_ReadCounts_deseq-Rlognormalized.txt", row.names=1)
exp_data <- read.delim("legend.txt", row.names=1)

## Set color scheme
mypalette<-brewer.pal(9,"Set1")

## Perform PCA, 2 possibilities shown, one for all samples, one for subset of samples. Both done on 10% most variant genes
## 1. PCA on all samples
raw_pca = prcomp(t(counts))

## Select 10% most variant genes
counts$variance <- apply(counts, MARGIN = 1, var)
variance <- sort(counts$variance, decreasing =T)
threshold <- variance[length(variance)*0.1] #10% most variant genes
counts_Top10V <- counts[which(counts$variance > threshold), 1:(ncol(counts)-1)]

raw_pca_Top10V = prcomp(t(counts_Top10V))

## Plotting PCA via ggbiplot
## Only PC1 versus PC2 possible
PCA_test_Top10V <- ggbiplot(raw_pca_Top10V, groups=exp_data$clone, obs.scale = 1, var.scale = 1 +
                              theme(axis.title.x=element_blank()))+
  theme(legend.direction = 'horizontal', legend.position = 'top')
PCA_test_Top10V

## 2. Do PCA for subset of samples
## Define subset, here d10 and d14
subset1 <- grep("day1",rownames(exp_data))
exp_data_subset <- exp_data[subset1,]
subset2 <- grep("day1",colnames(counts))
counts_subset <- counts[,subset2]

## Select 10% most variant genes
counts_subset$variance <- apply(counts_subset, MARGIN = 1, var)
variance_subset <- sort(counts_subset$variance, decreasing =T)
threshold_subset <- variance_subset[length(variance_subset)*0.1] #10% most variant genes
counts_Top10V_subset <- counts_subset[which(counts_subset$variance > threshold_subset), 1:(ncol(counts_subset)-1)]
raw_pca_Top10V_subset = prcomp(t(counts_Top10V_subset))

## Plotting sample subset PCA via ggbiplot
## Only PC1 versus PC2 possible
PCA_test_Top10V_subset <- ggbiplot(raw_pca_Top10V_subset, groups=exp_data_subset$genotype, obs.scale = 1, var.scale = 1 +
                              theme(axis.title.x=element_blank()))+
  theme(legend.direction = 'horizontal', legend.position = 'top')
PCA_test_Top10V_subset
