##---------------------------------------------------------------------------
## Author 	      : Markus Breunig (ORCID: 0000-0002-2817-1974) 
## Support        : Coding Templates from Thomas Engleitner (ORCID: 0000-0002-7989-681X)
## Date started 	: 2021
## Last modified  : 2025-03-24
## Aim 		        : Codes used for gene expression heatmaps from preprocessed and normalized RNA-seq bulk data
## Paper          : Breunig et al. 2025 - PLOS Genetics - Code generating Fig.5H, 6C-D, 7G
##--------------------------------------------------------------------------- 

## Load installed libraries
library("pheatmap")
library("biomaRt")
library("RColorBrewer")
library("tidyverse")
library("dplyr")

## Empty workspace
rm(list=ls())

## Load normalized expression data
Ndat <- read.table(file = "GSE282625_ReadCounts_deseq-Rlognormalized.txt", header=T, sep="\t", row.names=1)
targets = read.table("legend.txt",header=T,sep="\t", row.names=1)

##Duplicate matrix to overlay gene names to check presence of specific genes
Ndat2 = Ndat
Ndat2$GeneName <- rownames(Ndat)

## Set Heatmap color scheme
HeatColor = c("#053061","#0C3D74","#134B87","#1A599A","#2167AC","#2A72B2","#337EB8",
              "#3C89BE","#4795C4","#5BA2CB","#6FAFD2","#83BCD9","#96C7DF","#A6CFE3",
              "#B7D7E8","#C7E0ED","#D4E6F0","#DEEBF2","#E8F0F4","#F2F4F6","#F7F3F0","#F9ECE4","#FAE5D8","#FCDDCB","#FBD2BC","#F9C4AA","#F7B799","#F4A987",
              "#EE9878","#E6866A","#DF755D","#D7634F","#CE5146","#C53E3D","#BC2C34",
              "#B2192B","#A01228","#8D0C25","#7A0622","#67001F")

###--------------------------------------------------
## Different gene list sources are now described, all ending with a Genes vector comprised of the gene list of interest to be plotted
## 1. Blot Go-Term data sets
## Retrieve all HUGO gene symbols from investigated GO terms e.g.: 
##Execution phase of apoptosis: GO:0097194
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#go=c("GO:0097194")
#translation<- getBM(attributes= "hgnc_symbol",
#                    filters=c("go"),
#                    values=list(go), mart=human)
#colnames(translation) <- c("GeneName")
#overlap <- merge(translation,Ndat2)
#Genes = overlap$GeneName

## 2. Blot Reactome data sets
## Retrieve all HUGO gene symbols from reactome identifiers e.g.:
##Devolopmental Biology R-HSA-1266738
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#Reactome_ID=c("R-HSA-1266738")
#translation<- getBM(attributes= "hgnc_symbol",
#                    filters=c("reactome_gene"),
#                    values=list(Reactome_ID), mart=human)
#colnames(translation) <- c("GeneName")
#overlap <- merge(translation,Ndat2)
#Genes = overlap$GeneName
#translation <- unlist (translation, recursive =F) ##required without merge function to create vector input from list
#is.vector(translation)
#Genes = translation

## 3. Load and plot elsewhere saved gene/protein list
#Genes  <- read.delim2("ThisFileNameNeedsToBeAdded.txt", stringsAsFactors = F)

## Check if conversion to human gene nomenclature necessary, header, sep, data type might need to be adjusted
##Human to mouse conversion from following online workflow
#mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
#convert_mouse_to_human <- function(gene_list){
 # output = c()
 #for(gene in gene_list){
  #class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
 #if(!identical(class_key, integer(0)) ){
  #human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
 #for(human_gene in human_genes){
  #output = append(output,human_gene) }}}
#return (output)}

#Genes <- convert_mouse_to_human(Genes$GeneName) ##Column name might need to be adapted 
#Genes <- intersect(Genes,Ndat2$GeneName)

## 4. Manual Gene selection
Genes=c("PNPO", "MGMT", "CHCHD2", "TXNRD2", "TP53")
###--------------------------------------------------

## Actual plotting. Two possibilities shown
## 1. Plotting of Genes of interest without sample selection
#pheatmap(Ndat[Genes,],scale="row",cluster_cols=F,cluster_rows=F,col=HeatColor,clustering_method="ward.D2")

## 2. Plotting of Genes of interest wit sample selection and specification of plotting order
## Selection of covariate and samples 
Kovariate = "day" #"genotype" / "day"
Level = c( "0", "4", "10", "14") #"WT", "KO" / "0", "4", "10", "14"
Selection = rownames(targets[targets[,Kovariate] %in% Level,])
#order of display, please also specify then plot levels(Selection) instead of Selection in pheatmap function
##all samples
Selection <- factor(Selection, levels = c("X598_1_HUES8_Control_day0","X598_2_HUES8_Control_day0","X598_85_HUES8_Control_day0","X598_86_HUES8_Control_day0","X598_3_HUES8_EXDPF_KO_day0","X598_4_HUES8_EXDPF_KO_day0","X598_87_HUES8_EXDPF_KO_day0","X598_88_HUES8_EXDPF_KO_day0","X598_5_HUES8_Control_day4","X598_6_HUES8_Control_day4","X598_25_HUES8_Control_day4","X598_26_HUES8_Control_day4","X598_45_HUES8_Control_day4","X598_65_HUES8_Control_day4","X598_66_HUES8_Control_day4","X598_7_HUES8_EXDPF_KO_day4","X598_8_HUES8_EXDPF_KO_day4","X598_27_HUES8_EXDPF_KO_day4","X598_28_HUES8_EXDPF_KO_day4","X598_47_HUES8_EXDPF_KO_day4","X598_48_HUES8_EXDPF_KO_day4","X598_67_HUES8_EXDPF_KO_day4","X598_68_HUES8_EXDPF_KO_day4","X598_9_HUES8_Control_day10","X598_10_HUES8_Control_day10","X598_29_HUES8_Control_day10","X598_30_HUES8_Control_day10","X598_49_HUES8_Control_day10","X598_50_HUES8_Control_day10","X598_69_HUES8_Control_day10","X598_70_HUES8_Control_day10","X598_11_HUES8_EXDPF_KO_day10","X598_12_HUES8_EXDPF_KO_day10","X598_31_HUES8_EXDPF_KO_day10","X598_32_HUES8_EXDPF_KO_day10","X598_51_HUES8_EXDPF_KO_day10","X598_52_HUES8_EXDPF_KO_day10","X598_71_HUES8_EXDPF_KO_day10","X598_72_HUES8_EXDPF_KO_day10","X598_13_HUES8_Control_day14","X598_14_HUES8_Control_day14","X598_33_HUES8_Control_day14","X598_34_HUES8_Control_day14","X598_53_HUES8_Control_day14","X598_54_HUES8_Control_day14","X598_73_HUES8_Control_day14","X598_74_HUES8_Control_day14","X598_15_HUES8_EXDPF_KO_day14","X598_16_HUES8_EXDPF_KO_day14","X598_35_HUES8_EXDPF_KO_day14","X598_36_HUES8_EXDPF_KO_day14","X598_55_HUES8_EXDPF_KO_day14","X598_56_HUES8_EXDPF_KO_day14","X598_75_HUES8_EXDPF_KO_day14","X598_76_HUES8_EXDPF_KO_day14"))

## Plot Heatmap
pheatmap(Ndat[Genes,levels(Selection)],scale="row",cluster_cols=F,cluster_rows=T,col=HeatColor,clustering_method="ward.D2")

## Troubleshooting. Check complete cases, if clustering row does not work
# ForPlot <-  Ndat[Genes,] 
# ForPlot <- ForPlot[complete.cases(ForPlot),]
# pheatmap(ForPlot[,levels(Selection)],scale="row",cluster_cols=F,cluster_rows=T,col=HeatColor,clustering_method="ward.D2")

