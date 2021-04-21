#### calculates biadjacency matrix from vector of peptide sequences and vector
#### of corresponding protein accessions (separated by "/")

generate_01_matrix <- function(peptides, proteins) {
  require(pbapply)
  
  if(is.factor(peptides)) peptides <- as.character(peptides)
  if(is.factor(proteins)) proteins <- as.character(proteins)
  
  proteins_split <- strsplit(proteins, "/")
  
  protein_list <- sort(unique(unlist(proteins_split)))
  # Anzahl Proteine
  
  ### RES: biadjacency matrix
  ### peptides in rows
  ### proteins in columns
  ### TRUE -> peptide belongs to protein
  RES <- pbapply::pblapply(proteins_split, function(x) protein_list %in% x)
  RES <- t(as.data.frame(RES))
  colnames(RES) <- protein_list
  rownames(RES) <- peptides
  return(RES)
}


#### calculate connected components (submatrices) from biadjacency matrix
## M: biadjacency matrix
generate_submatrices <- function(M) {
  require(igraph)
  require(pbapply)
  
  G <- igraph::graph_from_incidence_matrix(matrix_01)
  Subgraphs <- igraph::decompose(G)
  Submatrix <- pbapply::pblapply(Subgraphs, as_incidence_matrix)
  return(list(Submatrix = Submatrix, Subgraphs = Subgraphs))
}



#### collapse protein nodes (merges duplicated columns of biadjacency submatrices)
## S: Submatrix list
## sparse: TRUE, if sparse matrix type is used
## fast: TRUE, if combining of column names should be skipped (makes calculation faster if there are large subgraphs)
form_proteingroups <- function(S, sparse = FALSE, fast = FALSE) {
  require(BBmisc)
  
  
  if (sparse) source("R Scripts/general_helper_functions/duplicated_for_sparse_matrices.R")
  
  for (i in seq_along(S)) { # for each submatrix in S
    print(i)
    tmp <- S[[i]]
    if(ncol(tmp) == 1) next
    
    if (sparse) {
      ind <- duplicated.dgCMatrix(tmp, MARGIN = 2)
      if (!any(ind)) next                                
      tmp2 <- tmp[, !ind, drop = FALSE]                    
      tmp3 <- tmp[, ind, drop = FALSE]
 
      if (!fast) {
        for (j in 1:ncol(tmp2)) {
          for (k in 1:ncol(tmp3)) {
            
            if (all(tmp3[,k] == tmp2[,j]))  {
              groupname <- BBmisc::collapse(c(colnames(tmp3)[k], colnames(tmp2)[j]), sep = ";")
              colnames(tmp2)[j] <- groupname
            }
          }
        }
      }
      
    } else {
      ind <- duplicated(tmp, MARGIN = 2)
      tmp2 <- tmp[, !ind, drop = FALSE]                     
      
      if (!fast) {
        for (j in 1:ncol(tmp2)) {
          ind <- apply(tmp, 2, function(x) all(x == tmp2[,j]))
          groupname <- BBmisc::collapse(colnames(tmp)[ind], sep = ";")
          colnames(tmp2)[j] <- groupname
        }
      }
      
    }
    S[[i]] <- tmp2
  }
  return(S)
}


#### adds peptide ratios to the list of submatrices,
## S: Submatrix
## fc: vector of fold changes
## peptides: list of peptides (to fc)
add_fc_to_submatrix <- function(S, fc, peptides) {
  S_new <- list()
  for (i in 1:length(S)) {
    S_i <- S[[i]]
    peptides_S <- rownames(S_i)
    ind <- match(peptides_S, peptides) 
    fc_S <- fc[ind]
    S_new[[i]] <- list(X = S_i, fc = fc_S)
  }
  return(S_new)
}



#### collapse peptide nodes
## S: submatrix list
## sparse: TRUE, if sparse matrix type is used
## fast: TRUE, if combining ofrow names should be skipped (makes calculation faster if there are large subgraphs)
## fc: TRUE, if S contains peptide ratios
merge_Peptides <- function(S, sparse = FALSE, fc = TRUE, fast = FALSE) {
  require(BBmisc)
  
  if (sparse) source("R Scripts/general_helper_functions/duplicated_for_sparse_matrices.R")
  
  if (fc) {
    for (i in seq_along(S)) {
      
      tmp <- S[[i]]$X
      fc <- S[[i]]$fc
      ind <- duplicated(tmp, MARGIN = 1)                   
      tmp2 <- tmp[!ind, , drop = FALSE]                      
      
      fc_tmp <- rep(NA, nrow(tmp2))
      for (j in 1:nrow(tmp2)) {
        ind <- apply(tmp, 1, function(x) all(x == tmp2[j,]))
        fc_tmp[j] <- 2^mean(log2(fc[ind]))
        groupname <- BBmisc::collapse(rownames(tmp)[ind], sep = ";")
        rownames(tmp2)[j] <- groupname
      }
      S[[i]]$X <- tmp2
      S[[i]]$fc <- fc_tmp
    }
  } else {
    for (i in seq_along(S)) {
      print(i)
      tmp <- S[[i]]
      if(sparse) {
        ind <- duplicated.dgCMatrix(tmp, MARGIN = 1)
        
        if (!any(ind)) next                                   
        
        tmp2 <- tmp[!ind, , drop = FALSE]       
        tmp3 <- tmp[ind, , drop = FALSE]
        
        if (!fast) {
          for (j in 1:nrow(tmp2)) {
            for (k in 1:nrow(tmp3)) {
              
              if (all(tmp3[k,] == tmp2[j,]))  {
                groupname <- BBmisc::collapse(c(rownames(tmp3)[k], rownames(tmp2)[j]), sep = ";")   
                rownames(tmp2)[j] <- groupname
              }
            }
          }
        }
        
        S[[i]]<- tmp2
        
      } else {
        ind <- duplicated(tmp, MARGIN = 1)  
        tmp2 <- tmp[!ind, , drop = FALSE]                    
        
        fc_tmp <- rep(NA, nrow(tmp2))
        for (j in 1:nrow(tmp2)) {
          ind <- apply(tmp, 1, function(x) all(x == tmp2[j,]))
          groupname <- BBmisc::collapse(rownames(tmp)[ind], sep = ";")  
          rownames(tmp2)[j] <- groupname
        }
        S[[i]] <- tmp2
      }
    }
  }
  return(S)
}


