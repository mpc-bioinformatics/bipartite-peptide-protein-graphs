library(BBmisc)    # for collapse()
library(seqinr)    # for reading in fasta files
library(limma)     # for strsplit2()
library(pbapply)   # for progress bars for apply functions

source("helper_function/assign_protein_accessions.R")
source("helper_function/Digest2.R")
source("helper_function/calculate_biadjacency_matrix_and_submatrices.R")
source("helper_function/calculate_peptide_ratios.R")

DATA <- read.table("data/D2/D2_quant/peptides.txt", 
                    sep = "\t", 
                    header = TRUE, stringsAsFactors = FALSE)

####  extract peptide intensities and rename columns ####
LFQ_values <- DATA[, grepl("LFQ", colnames(DATA))]
# order by concentration of spiked-in UPS1 standard:
LFQ_values <- LFQ_values[, c(25:27,  4:6,   13:15,  22:24, 10:12, 19:21, 1:3,  7:9, 16:18)]
colnames(LFQ_values) <- paste0("state", rep(1:9, each = 3), "_", rep(1:3, 9))


################################################################################
#### Step 1) remove decoy proteins ####

LFQ_values <- LFQ_values[DATA$Reverse == "", ]
DATA <- DATA[DATA$Reverse == "", ] 

################################################################################
#### Step 2) summarize modifications (given in round brackets) ####

################################################################################
#### Step 3) assign peptides to protein accessions via in silico digestion of the fasta file ####

sequence <- DATA$Sequence
fasta <- read.fasta(file = "data/D2/D2_fasta/2020_01_31_proteome_S_cerevisae_with_isoforms.fasta", 
                    seqtype = "AA", as.string = TRUE)
fasta_vec <- unlist(fasta)
protein_accessions <- strsplit2(attr(fasta, "name"), "\\|")[,2]
fasta_UPS1 <- read.fasta(file = "data/D2/D2_fasta/ups1-ups2-sequences.fasta", 
                         seqtype = "AA", as.string = TRUE)
fasta_vec <- c(fasta_vec, unlist(fasta_UPS1))
protein_accessions <- c(protein_accessions, strsplit2(attr(fasta_UPS1, "name"), "\\|")[,1])
fasta_contaminant <- read.fasta(file = "data/D2/D2_fasta/contaminants.fasta", 
                                seqtype = "AA", as.string = TRUE)
fasta_vec <- c(fasta_vec, unlist(fasta_contaminant))
protein_accessions <- c(protein_accessions, paste0("CON_", attr(fasta_contaminant, "name")))
names(fasta_vec) <- protein_accessions

peptides <- DATA$Sequence
proteins <- assign_protein_accessions(peptides, fasta_vec)

D <- cbind(peptide = peptides, protein = proteins, LFQ_values)
write.table(D, "data/D2/D2_quant/preprocessed_peptide_data_D2.txt", 
            sep = "\t", row.names = FALSE)

################################################################################
#### Step 4) calculate 0/1 biadjacency matrix matrix ####

DATA <- read.table(file = "data/D2/D2_quant/preprocessed_peptide_data_D2.txt", sep = "\t", header = TRUE)
matrix_01 <- generate_01_matrix(DATA$peptide, DATA$protein)
saveRDS(matrix_01, "data/D2/D2_quant/preprocessed/01Matrix.rds")


################################################################################
#### Step 5) find sub matrices / connected components ####

matrix_01 <- readRDS("data/D2/D2_quant//preprocessed/01Matrix.rds")
RES <- generate_submatrices(matrix_01)
saveRDS(RES$Submatrix, file = "data/D2/D2_quant//preprocessed/Submatrix.rds")
saveRDS(RES$Subgraphs, file = "data/D2/D2_quant//preprocessed/Subgraphs.rds")


################################################################################
#### Step 6) collapse protein nodes ####

S <- readRDS("data/D2/D2_quant/preprocessed/Submatrix.rds")
S_proteingroups <- form_proteingroups(S$Submatrix)
saveRDS(S_proteingroups, file = "data/D2/D2_quant/preprocessed/Submatrix_proteinGroups.rds")


################################################################################
#### Step 7) calculate peptide ratios / fold changes ####

D <- read.table(file = "data/D2/D2_quant/preprocessed_peptide_data_D2.txt", header = TRUE)
D_mean_NA_mind2 <- aggregate_replicates(D, method = "mean", use0 = FALSE, missing.limit = 0.4, 
                                        group = factor(substr(colnames(D)[-(1:2)], 1, 6)), accession.cols = 1:2)
write.table(D_mean_NA_mind2, file = "data/D2/D2_quant/preprocessed/D_mean_NA_mind2.txt", 
            row.names = FALSE, sep = "\t")

### one missing value allowed:
D_aggr <- read.table("data/D2/D2_quant/preprocessed/D_mean_NA_mind2.txt", sep = "\t", header = TRUE)
FC_1_2 <- foldChange(D_aggr, "state1", "state2")
FC_1_3 <- foldChange(D_aggr, "state1", "state3")
FC_1_4 <- foldChange(D_aggr, "state1", "state4")
FC_1_5 <- foldChange(D_aggr, "state1", "state5")
FC_1_6 <- foldChange(D_aggr, "state1", "state6")
FC_1_7 <- foldChange(D_aggr, "state1", "state7")
FC_1_8 <- foldChange(D_aggr, "state1", "state8")
FC_1_9 <- foldChange(D_aggr, "state1", "state9")
FC_2_3 <- foldChange(D_aggr, "state2", "state3")
FC_2_4 <- foldChange(D_aggr, "state2", "state4")
FC_2_5 <- foldChange(D_aggr, "state2", "state5")
FC_2_6 <- foldChange(D_aggr, "state2", "state6")
FC_2_7 <- foldChange(D_aggr, "state2", "state7")
FC_2_8 <- foldChange(D_aggr, "state2", "state8")
FC_2_9 <- foldChange(D_aggr, "state2", "state9")
FC_3_4 <- foldChange(D_aggr, "state3", "state4")
FC_3_5 <- foldChange(D_aggr, "state3", "state5")
FC_3_6 <- foldChange(D_aggr, "state3", "state6")
FC_3_7 <- foldChange(D_aggr, "state3", "state7")
FC_3_8 <- foldChange(D_aggr, "state3", "state8")
FC_3_9 <- foldChange(D_aggr, "state3", "state9")
FC_4_5 <- foldChange(D_aggr, "state4", "state5")
FC_4_6 <- foldChange(D_aggr, "state4", "state6")
FC_4_7 <- foldChange(D_aggr, "state4", "state7")
FC_4_8 <- foldChange(D_aggr, "state4", "state8")
FC_4_9 <- foldChange(D_aggr, "state4", "state9")
FC_5_6 <- foldChange(D_aggr, "state5", "state6")
FC_5_7 <- foldChange(D_aggr, "state5", "state7")
FC_5_8 <- foldChange(D_aggr, "state5", "state8")
FC_5_9 <- foldChange(D_aggr, "state5", "state9")
FC_6_7 <- foldChange(D_aggr, "state6", "state7")
FC_6_8 <- foldChange(D_aggr, "state6", "state8")
FC_6_9 <- foldChange(D_aggr, "state6", "state9")
FC_7_8 <- foldChange(D_aggr, "state7", "state8")
FC_7_9 <- foldChange(D_aggr, "state7", "state9")
FC_8_9 <- foldChange(D_aggr, "state8", "state9")

FC_list <- mget(ls(pattern = "FC"))

FC <- do.call(cbind, FC_list)
FC <- as.data.frame(FC)
FC <- cbind(D_aggr[, 1:2], FC)
write.table(FC, file = "data/D2/D2_quant/preprocessed/FC_NA_mind2.txt", row.names = FALSE, sep = "\t")


################################################################################
#### Step 8) re-calculate sub matrices but only with peptides with valid ratios per pairwise comparison #### 

FC.table <- read.table(file = "data/D2/D2_quant/preprocessed/FC_NA_mind2.txt", header = TRUE, sep = "\t")

system.time({
for (i in 1:8) {
  for (j in (i + 1):9) {
    print(c(i,j))
    
    fc <- FC.table[, paste0("FC_", i, "_", j)]
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
  }
}
})

save(list =  ls(pattern = "Submatrix"), 
     file = "data/D2/D2_quant/preprocessed/Submatrices.RData")



################################################################################
#### Step 9) collapse peptide nodes ####

load("data/D2/D2_quant/preprocessed/Submatrices.RData")

system.time({
for (i in 1:8) {
  for (j in (i + 1):9) {
    S <- get(paste("Submatrix", i, j, sep = "_"))
    S_new <- merge_Peptides(S)
    assign(paste0("Submatrix_merged_Peptides_", i, "_", j), S_new)
  }
}
})
  
save(list =  ls(pattern = "Submatrix_merged"), 
     file = "data/D2/D2_quant/preprocessed/Submatrices_merged_Peptides.RData")

