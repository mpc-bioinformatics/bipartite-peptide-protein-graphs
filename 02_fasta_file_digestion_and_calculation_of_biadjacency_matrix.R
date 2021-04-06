library(seqinr)    # for reading in fasta files
library(pbapply)   # for progress bars for apply functions
library(Matrix)    # for sparse matrices
library(limma)     # for strsplit2()

source("helper_functions/Digest2.R")


################################################################################
#### read in fasta files and compute complete biadjacency matrix ####

################################################################################
### D1_fasta

fasta <- read.fasta(file = "data/D1/D1_fasta/uniprot-proteome-mus_musculus-spiked_accessions-cRAP-iRT-2017_12.fasta", 
                    seqtype = "AA", as.string = TRUE)

digested_proteins <- pblapply(fasta, function(x) {
  sequ <- x
  class(sequ) <- NULL
  y <- try({Digest2(sequ, enzyme = "trypsin", missed = 2)})
  ind <- nchar(as.character(y$sequence)) >= 5 & nchar(as.character(y$sequence)) <= 50 # at least 5 AA and at most 50 AA
  return(as.character(y$sequence[ind]))  
})

peptides <- sort(unique(unlist(digested_proteins)))  # character of peptide sequences, ordered alphabetically
peptides_fac <- factor(peptides)

#### construction of biadjacency matrix as a sparse matrix
sparseM <- Matrix(0, nrow = length(peptides), ncol = length(digested_proteins), sparse = TRUE)
for (i in 1:length(digested_proteins)) { 
  print(i)
  ind <- which(peptides %in% digested_proteins[[i]])
  if (length(ind) > 0) sparseM[ind, i] <- 1
}

colnames(sparseM) <- names(digested_proteins)
rownames(sparseM) <- peptides

save(sparseM, file = "data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.RData")
saveRDS(sparseM, file = "data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.rds")


################################################################################
### D2_fasta

fasta1 <- read.fasta(file = "data/D2/D2_fasta/2020_01_31_proteome_S_cerevisae_with_isoforms.fasta", 
                    seqtype = "AA", as.string = TRUE)
names(fasta1) <- strsplit2(names(fasta1), "\\|")[,2]
fasta2 <- read.fasta(file = "data/D2/D2_fasta/ups1-ups2-sequences.fasta", 
                     seqtype = "AA", as.string = TRUE)
names(fasta2) <- strsplit2(names(fasta2), "\\|")[,1]
fasta3 <- read.fasta(file = "data/D2/D2_fasta/MaxQuant_contaminants_downloaded_20200527.fasta", 
                     seqtype = "AA", as.string = TRUE)
fasta <- c(fasta1, fasta2, fasta3)


digested_proteins <- pblapply(fasta, function(x) {
  sequ <- x
  class(sequ) <- NULL
  y <- try({Digest2(sequ, enzyme = "trypsin", missed = 2)})
  ind <- nchar(as.character(y$sequence)) >= 5 & nchar(as.character(y$sequence)) <= 50 # at least 5 AA and at most 50 AA
  return(as.character(y$sequence[ind]))  
})


peptides <- sort(unique(unlist(digested_proteins)))
peptides_fac <- factor(peptides)

sparseM <- Matrix(0, nrow = length(peptides), ncol = length(digested_proteins), sparse = TRUE)
for (i in 1:length(digested_proteins)) { ### 1
  print(i)
  ind <- which(peptides %in% digested_proteins[[i]])
  if (length(ind) > 0) sparseM[ind, i] <- 1
}

colnames(sparseM) <- names(digested_proteins)
rownames(sparseM) <- peptides

save(sparseM, file = "data/D2/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.RData")
saveRDS(sparseM, file = "data/D2/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.rds")





