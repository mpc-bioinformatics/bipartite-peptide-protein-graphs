### installation of required packages

if(!require(seqinr)) install.packages("seqinr")
if(!require(OrgMassSpecR)) install.packages("OrgMassSpecR")
if(!require(igraph)) install.packages("igraph")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(reshape2)) install.packages("reshape2")
if(!require(BBmisc)) install.packages("BBmisc")
if(!require(pbapply)) install.packages("pbapply")
if(!require(Matrix)) install.packages("Matrix")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(xtable)) install.packages("xtable")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(cowplot)) install.packages("cowplot")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(limma)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
}



################################################################################
### create necessary folder structure to store (intermediate) results

dir.create("data/D1/D1_fasta/preprocessed/", recursive = TRUE)
dir.create("data/D1/D1_fasta/isomorph_classes/", recursive = TRUE)
dir.create("data/D1/D1_quant/preprocessed/", recursive = TRUE)
dir.create("data/D1/D1_quant/isomorph_classes/", recursive = TRUE)

dir.create("data/D2_without_isoforms/D2_fasta/preprocessed/", recursive = TRUE)
dir.create("data/D2_without_isoforms/D2_fasta/isomorph_classes/", recursive = TRUE)
dir.create("data/D2_without_isoforms/D2_quant/preprocessed/", recursive = TRUE)
dir.create("data/D2_without_isoforms/D2_quant/isomorph_classes/", recursive = TRUE)

dir.create("data/D3_without_isoforms/D3_fasta/preprocessed/", recursive = TRUE)
dir.create("data/D3_without_isoforms/D3_fasta/isomorph_classes/", recursive = TRUE)
dir.create("data/D3_without_isoforms/D3_quant/preprocessed/", recursive = TRUE)
dir.create("data/D3_without_isoforms/D3_quant/isomorph_classes/", recursive = TRUE)

dir.create("data/D3/D3_fasta/preprocessed/", recursive = TRUE)
dir.create("data/D3/D3_fasta/isomorph_classes/", recursive = TRUE)
dir.create("data/D3/D3_quant/preprocessed/", recursive = TRUE)
dir.create("data/D3/D3_quant/isomorph_classes/", recursive = TRUE)

dir.create("Paper/Paper 1/figures/supplement/", recursive = TRUE)


