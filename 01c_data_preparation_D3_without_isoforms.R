library(BBmisc)    # for collapse()
library(seqinr)    # for reading in fasta files
library(limma)     # for strsplit2()
library(pbapply)   # for progress bars for apply functions
library(vsn)      # lts-normalization
library(tidyverse)


source("helper_functions/assign_protein_accessions.R")
source("helper_functions/Digest2.R")
source("helper_functions/calculate_biadjacency_matrix_and_submatrices.R")
source("helper_functions/calculate_peptide_ratios.R")

################################################################################

DATA <- read.table("data/D3_without_isoforms/D3_quant/peptides.txt",
                   sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE)

####  extract peptide intensities and rename columns ####
LFQ_values <- DATA[, grepl("LFQ", colnames(DATA))]
LFQ_values <- LFQ_values[, c(4:6, 1:3)]
colnames(LFQ_values) <- paste0("state", rep(1:2, each = 3), "_", rep(1:3, 2))


################################################################################
#### Step 1) remove decoy proteins ####

ind_notreverse <- DATA$Reverse == ""

LFQ_values <- LFQ_values[ind_notreverse, ]
LFQ_values[LFQ_values == 0] <- NA

DATA <- DATA[ind_notreverse,]


################################################################################
#### Step 2) summarize modifications (given in round brackets) ####

################################################################################
#### Step 3) assign peptides to protein accessions via in silico digestion of the fasta file ####

sequence <- DATA$Sequence
fasta <- read.fasta(file = "data/D3_without_isoforms/D3_fasta/20210413_uniprot-proteome UP000005640_HomoSapiens_v2021_02.fasta",
                    seqtype = "AA", as.string = TRUE)
fasta_vec <- unlist(fasta)
protein_accessions <- strsplit2(attr(fasta, "name"), "\\|")[,2]

fasta_EColi <- read.fasta(file = "data/D3_without_isoforms/D3_fasta/20210413_uniprot-proteome UP000000625_EColi_v2021_02.fasta",
                         seqtype = "AA", as.string = TRUE)
fasta_vec <- c(fasta_vec, unlist(fasta_EColi))
protein_accessions <- c(protein_accessions, strsplit2(attr(fasta_EColi, "name"), "\\|")[,2])
fasta_contaminant <- read.fasta(file = "data/D3_without_isoforms/D3_fasta/MaxQuant_contaminants_downloaded_20200527.fasta",
                                seqtype = "AA", as.string = TRUE)
fasta_vec <- c(fasta_vec, unlist(fasta_contaminant))
protein_accessions <- c(protein_accessions, paste0("CON_", attr(fasta_contaminant, "name")))
names(fasta_vec) <- protein_accessions

peptides <- DATA$Sequence
proteins <- assign_protein_accessions(peptides, fasta_vec)


D <- cbind(peptide = peptides, protein = proteins, LFQ_values)
D <- D[nchar(D$peptide) >= 7 & nchar(D$peptide) <= 50, ]
write.table(D, "data/D3_without_isoforms/D3_quant/preprocessed_peptide_data_D3_without_isoforms.txt",
            sep = "\t", row.names = FALSE)


################################################################################
#### Step 4) calculate 0/1 biadjacency matrix matrix ####

DATA <- read.table(file = "data/D3_without_isoforms/D3_quant/preprocessed_peptide_data_D3_without_isoforms.txt", sep = "\t", header = TRUE)
matrix_01 <- generate_01_matrix(DATA$peptide, DATA$protein)  ### approx. 1 minute
saveRDS(matrix_01, "data/D3_without_isoforms/D3_quant/preprocessed/01Matrix_without_isoforms.rds")


################################################################################
#### Step 5) find sub matrices / connected components ####

matrix_01 <- readRDS("data/D3_without_isoforms/D3_quant/preprocessed/01Matrix_without_isoforms.rds")
RES <- generate_submatrices(matrix_01)
saveRDS(RES$Submatrix, file = "data/D3_without_isoforms/D3_quant/preprocessed/Submatrix_without_isoforms.rds")
saveRDS(RES$Subgraphs, file = "data/D3_without_isoforms/D3_quant/preprocessed/Subgraphs_without_isoforms.rds")


################################################################################
#### Step 6) collapse protein nodes ####

S <- readRDS("data/D3_without_isoforms/D3_quant/preprocessed/Submatrix_without_isoforms.rds")
S_proteingroups <- form_proteingroups(S)
saveRDS(S_proteingroups, file = "data/D3_without_isoforms/D3_quant/preprocessed/Submatrix_proteinGroups_without_isoforms.rds")


################################################################################
#### Step 7) calculate peptide ratios / fold changes ####

D <- read.table(file = "data/D3_without_isoforms/D3_quant/preprocessed_peptide_data_D3_without_isoforms.txt", header = TRUE)
D_mean_NA_mind2 <- aggregate_replicates(D, method = "mean", use0 = FALSE, missing.limit = 0.4,
                                        group = factor(substr(colnames(D)[-c(1:2, 9)], 1, 6)), accession.cols = c(1:2,9))
write.table(D_mean_NA_mind2, file = "data/D3_without_isoforms/D3_quant/preprocessed/D_mean_NA_mind2.txt",
            row.names = FALSE, sep = "\t")

### one missing value allowed:
D_aggr <- read.table("data/D3_without_isoforms/D3_quant/preprocessed/D_mean_NA_mind2.txt", sep = "\t", header = TRUE)
FC_1_2 <- foldChange(D_aggr, "state1", "state2")

FC <- FC_1_2
FC <- as.data.frame(FC)
FC <- cbind(D_aggr[, 1:2], FC)
write.table(FC, file = "data/D3_without_isoforms/D3_quant/preprocessed/FC_NA_mind2.txt", row.names = FALSE, sep = "\t")


################################################################################
#### Step 8) re-calculate sub matrices but only with peptides with valid ratios per pairwise comparison ####


FC.table <- read.table(file = "data/D3_without_isoforms/D3_quant/preprocessed/FC_NA_mind2.txt", header = TRUE, sep = "\t")

i = 1
j = 2

print(c(i,j))

fc <- FC.table$FC
peptide <- as.character(FC.table$peptide)
peptide <- peptide[!is.na(fc)]
protein <- as.character(FC.table$protein)
protein <- protein[!is.na(fc)]
fc <- na.omit(fc)

matrix_01 <- generate_01_matrix(peptide, protein)

submat <- generate_submatrices(matrix_01)$Submatrix
submat <- form_proteingroups(submat)

submat_fc <- add_fc_to_submatrix(submat, fc, peptide)

assign(paste0("Submatrix_", i, "_", j), submat_fc)


save(list =  ls(pattern = "Submatrix"),
     file = "data/D3_without_isoforms/D3_quant/preprocessed/Submatrices.RData")



################################################################################
#### Step 9) collapse peptide nodes ####

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices.RData")

system.time({
  for (i in 1) {
    for (j in 2) {
      S <- get(paste("Submatrix", i, j, sep = "_"))
      S_new <- merge_Peptides(S)
      assign(paste0("Submatrix_merged_Peptides_", i, "_", j), S_new)
    }
  }
})

save(list =  ls(pattern = "Submatrix_merged"),
     file = "data/D3_without_isoforms/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")

