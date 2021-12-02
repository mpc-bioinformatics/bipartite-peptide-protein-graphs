library(BBmisc)   # for collapse()
library(seqinr)   # for reading in fasta files
library(pbapply)  # for progress bars for apply functions

source("helper_functions/assign_protein_accessions.R")
source("helper_functions/Digest2.R")
source("helper_functions/generate_01_matrix_and_submatrices.R")
source("helper_functions/calculate_peptide_ratios.R")

DATA <- read.table("data/D1/D1_quant/quantified_peptides-featureFinderCentroided.csv",
           sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 3)

Peptides <- DATA$peptide
Proteins <- DATA$protein


################################################################################
#### Step 1) remove decoy proteins ("s" before the accession) ####

Proteins_split <- strsplit(Proteins, "/")
Proteins_split2 <- list()

for (i in 1:nrow(DATA)) {
  x <- Proteins_split[[i]]
  decoy_ind <- grep("s", x, ignore.case = FALSE)
  if(any(regexpr("s", x, ignore.case = FALSE) > 1)) { # only remove accession with "s" as first letter
    decoy_ind <- decoy_ind[decoy_ind != which(regexpr("s", x, ignore.case = FALSE) > 1)]
  }
  if(length(decoy_ind) == 0) {
    Proteins_split2[[i]] <- x
  } else {
    Proteins_split2[[i]] <- x[-decoy_ind]
  }
}

### remove peptides that do not belong to any proteins after decoy removal:
ind_decoyPeptide <- sapply(Proteins_split2, length)==0
DATA <- DATA[!ind_decoyPeptide, ]
Proteins_split2 <- Proteins_split2[!ind_decoyPeptide]
Peptides <- Peptides[!ind_decoyPeptide]
Proteins <- sapply(Proteins_split2, function(x) BBmisc::collapse(x, sep = "/"))
DATA$protein <- Proteins
DATA <- DATA[, -c(3,4,20)]

write.table(DATA, file = "data/D1/D1_quant/preprocessed/datatable_without_decoys.txt", sep = "\t", row.names = FALSE)


################################################################################
#### Step 2) summarize modifications (given in round brackets) ####

DATA <- read.table(file = "data/D1/D1_quant/preprocessed/datatable_without_decoys.txt", sep = "\t", header = TRUE)

Peptides_sequence <- gsub( " *\\(.*?\\) *", "", DATA$peptide) # removes brackets and everything between them
Peptides_sequence <- gsub(" ", "", Peptides_sequence, fixed = TRUE) # removes white space
Peptides_sequence <- gsub(".", "", Peptides_sequence, fixed = TRUE) # remoces dot

ind <- duplicated(Peptides_sequence)  # duplicate = has multiple modifies versions
proteins_ <- DATA$protein[!ind]

DATA_aggregated_modifications <- DATA
DATA_aggregated_modifications$peptide <- Peptides_sequence
DATA_aggregated_modifications[is.na(DATA_aggregated_modifications)] <- 0  # set NAs to zero

### sum up intensities of same peptide sequences (w/o modifications) + add column with protein accessions:
DATA_aggregated_modifications <- aggregate(. ~ peptide, DATA_aggregated_modifications[, -2],
                                           FUN = function(x){sum(x, na.rm = TRUE)}, drop = FALSE)
DATA_aggregated_modifications <- cbind(peptide = DATA_aggregated_modifications[,1],
                                       protein = proteins_[order(unique(Peptides_sequence))], DATA_aggregated_modifications[, -1])

write.table(DATA_aggregated_modifications, "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications.txt", row.names = FALSE, sep = "\t")


################################################################################
#### Step 3) assign peptides to protein accessions via in silico digestion of the fasta file ####

DATA <- read.table("data/D1/D1_quant/preprocessed/datatable_aggregated_modifications.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)

sequence <- DATA$peptide
fasta <- read.fasta(file = "data/D1/D1_fasta/uniprot-proteome-mus_musculus-spiked_accessions-cRAP-iRT-2017_12.fasta",
                    seqtype = "AA", as.string = TRUE)
fasta_vec <- unlist(fasta)
fasta_protein_acc <- attr(fasta, "name")

proteins <- assign_protein_accessions(sequence, fasta_vec)

DATA$protein <- proteins
write.table(DATA, "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications2.txt",
            row.names = FALSE, sep = "\t")


DATA <- read.table("data/D1/D1_quant/preprocessed/datatable_aggregated_modifications2.txt", sep = "\t", header = TRUE)
DATA <- DATA[DATA$protein != "", ]  # remove peptides that could not be assigned to a protein
DATA <- DATA[nchar(DATA$peptide) >= 7 & nchar(DATA$peptide) <= 50, ]  # remove peptides outside of the desired length (7-50 AA)
write.table(DATA, "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications3.txt",
            row.names = FALSE, sep = "\t")


################################################################################
#### Step 4) calculate 0/1 biadjacency matrix matrix ####

### all peptides that are assigned to at least one protein:
DATA <- read.table(file = "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications2.txt",
                   sep = "\t", header = TRUE)
DATA <- DATA[DATA$protein != "", ]
matrix_01 <- generate_01_matrix(DATA$peptide, DATA$protein)
saveRDS(matrix_01, "data/D1/D1_quant/preprocessed/01Matrix_all_valid_peptides.rds")

### only peptides between size of 7 and 50:
DATA <- read.table(file = "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications3.txt",
                   sep = "\t", header = TRUE)
matrix_01 <- generate_01_matrix(DATA$peptide, DATA$protein)
saveRDS(matrix_01, "data/D1/D1_quant/preprocessed/01Matrix.rds")


################################################################################
#### Step 5) find sub matrices / connected components ####

matrix_01 <- readRDS("data/D1/D1_quant/preprocessed/01Matrix.rds")
Submatrix <- generate_submatrices(matrix_01)
saveRDS(Submatrix, file = "data/D1/D1_quant/preprocessed/Submatrix.rds")


################################################################################
#### Step 6) collapse protein nodes ####

S <- readRDS("data/D1/D1_quant/preprocessed/Submatrix.rds")
S_proteingroups <- form_proteingroups(S$Submatrix)
saveRDS(S_proteingroups, file = "data/D1/D1_quant/preprocessed/Submatrix_proteinGroups.rds")


################################################################################
#### Step 7) calculate peptide ratios / fold changes ####

D <- read.table(file = "data/D1/D1_quant/preprocessed/datatable_aggregated_modifications3.txt", header = TRUE)
D_mean_NA_mind2 <- aggregate_replicates(D, method = "mean", use0 = FALSE, missing.limit = 0.4,
                                        group = factor(rep(paste0("state", 1:5), 3)))
write.table(D_mean_NA_mind2, file = "data/D1/D1_quant/preprocessed/D_mean_NA_mind2.txt", row.names = FALSE, sep = "\t")

D_aggr <- read.table("data/D1/D1_quant/preprocessed/D_mean_NA_mind2.txt", sep = "\t", header = TRUE)
FC_1_2 <- foldChange(D_aggr, "state1", "state2")
FC_1_3 <- foldChange(D_aggr, "state1", "state3")
FC_1_4 <- foldChange(D_aggr, "state1", "state4")
FC_1_5 <- foldChange(D_aggr, "state1", "state5")
FC_2_3 <- foldChange(D_aggr, "state2", "state3")
FC_2_4 <- foldChange(D_aggr, "state2", "state4")
FC_2_5 <- foldChange(D_aggr, "state2", "state5")
FC_3_4 <- foldChange(D_aggr, "state3", "state4")
FC_3_5 <- foldChange(D_aggr, "state3", "state5")
FC_4_5 <- foldChange(D_aggr, "state4", "state5")
FC <- cbind(FC_1_2, FC_1_3, FC_1_4, FC_1_5, FC_2_3, FC_2_4, FC_2_5, FC_3_4, FC_3_5, FC_4_5)
FC <- as.data.frame(FC)
FC <- cbind(D_aggr[, 1:2], FC)
write.table(FC, file = "data/D1/D1_quant/preprocessed/FC_NA_mind2.txt", row.names = FALSE, sep = "\t")


################################################################################
#### Step 8) re-calculate sub matrices but only with peptides with valid ratios per pairwise comparison ####

FC.table <- read.table(file = "data/D1/D1_quant/preprocessed/FC_NA_mind2.txt", header = TRUE, sep = "\t")

for (i in 1:4) {
  for (j in (i + 1):5) {
    print(c(i,j))
    fc <- FC.table[, paste0("FC_", i, "_", j)]
    peptide <- as.character(FC.table$peptide)
    peptide <- peptide[!is.na(fc)]
    protein <- as.character(FC.table$protein)
    protein <- protein[!is.na(fc)]
    fc <- na.omit(fc)
    matrix_01 <- generate_01_matrix(peptide, protein)
    submat <- generate_submatrices(matrix_01)
    submat <- form_proteingroups(submat)
    submat_fc <- add_fc_to_submatrix(submat, fc, peptide)
    assign(paste0("Submatrix_", i, "_", j), submat_fc)
  }
}


save(Submatrix_1_2,
     Submatrix_1_3,
     Submatrix_1_4,
     Submatrix_1_5,
     Submatrix_2_3,
     Submatrix_2_4,
     Submatrix_2_5,
     Submatrix_3_4,
     Submatrix_3_5,
     Submatrix_4_5,
     file = "data/D1/D1_quant/preprocessed/Submatrices.RData")

################################################################################
#### Step 9) collapse peptide nodes ####

load("data/D1/D1_quant//preprocessed/Submatrices.RData")
for (i in 1:4) {
  for (j in (i + 1):5) {
    S <- get(paste("Submatrix", i, j, sep = "_"))
    S_new <- merge_Peptides(S)
    assign(paste0("Submatrix_merged_Peptides_", i, "_", j), S_new)
  }
}

save(Submatrix_merged_Peptides_1_2,
  Submatrix_merged_Peptides_1_3,
  Submatrix_merged_Peptides_1_4,
  Submatrix_merged_Peptides_1_5,
  Submatrix_merged_Peptides_2_3,
  Submatrix_merged_Peptides_2_4,
  Submatrix_merged_Peptides_2_5,
  Submatrix_merged_Peptides_3_4,
  Submatrix_merged_Peptides_3_5,
  Submatrix_merged_Peptides_4_5,
  file = "data/D1/D1_quant/\\preprocessed\\Submatrices_merged_Peptides.RData")
