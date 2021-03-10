#### assign protein accessions to peptide sequences
## sequence: vector of peptide sequences
## fasta_vec: names vector of protein sequences from fasta file(s) 

assign_protein_accessions <- function(sequence, fasta_vec) {
  require(pbapply)
  require(BBmisc)
  source("R Scripts/Digest2.R")
  
  protein_accessions <- names(fasta_vec)
  
  assigned_proteins <- pblapply(sequence, function(x) {
    ind <- which(grepl(x, fasta_vec))
    proteins <- protein_accessions[ind]
    is_tryptic <- rep(TRUE, length(ind)) 
    for (i in 1:length(ind)) {
      dig <- Digest2(fasta_vec[ind[i]], missed = 2)
      if (!(x %in% dig$sequence)) {
        if(any(dig$start == 1)) {
          dig2 <- dig[dig$start == 1, ]  # only peptides at protein N-terminus, where initial M could possible have been cut off
          if (!(paste0("M", x) %in% dig2$sequence)) {
            is_tryptic[i] <- FALSE
            next()
          }
        } else {
          is_tryptic[i] <- FALSE
        }
      }
    }
    proteins <- proteins[is_tryptic]
    proteins <- BBmisc::collapse(proteins, sep = "/")
    return(proteins)
  })
  
  return(unlist(assigned_proteins))
}


