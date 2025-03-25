##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Date started 	: 2022-06-24
## Last modified  : 2025-03-24
## Aim 		        : Codes used for variance heatmap/hierarchical clustering from preprocessed and normalized RNA-seq bulk data
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.5C
##--------------------------------------------------------------------------- 

## Load installed libraries
library("RColorBrewer")
library("pheatmap")

## Empty workspace
rm(list=ls())

## Load data
counts <- read.delim("GSE282625_ReadCounts_deseq-Rlognormalized.txt", row.names=1)
legend <- read.delim("legend.txt", row.names=1)
legend$summary <- paste(legend$clone, sep = "_day", legend$day)

## A distance matrix gives us an overview over similarities and dissimilarities between samples. 
sampleDists <- dist(t(counts))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(legend$summary)
colnames(sampleDistMatrix) <- paste(legend$Differentiation)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## Plotting of heatmap. We have to provide a hierarchical clustering to the heatmap function based on the sample distances, 
## or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.
## default is Euclidean distance
pheatmap(sampleDistMatrix, cutree_rows = 5, cutree_cols = 5,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

