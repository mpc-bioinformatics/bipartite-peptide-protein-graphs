#### aggregates replicates of the same experimental group
# X: data set
# missing.limit:  proportion of missing values that is allowed (e.g. 0 means no missings allowed)
# method = "mean", "sum" oder "median"
# use0: TRUE: calculate with 0, FALSE: replace 0 with NA
# group: groups for aggregation as a factor
# accessions.cols: columns with accession etc that are skipped for aggregation
aggregate_replicates <- function(X, missing.limit = 0, method = "mean", use0 = FALSE, 
                                 group = factor(substr(colnames(D)[-(1:2)], 1, 6)), accession.cols = 1:2) {
  require(robustbase)
  
  acc <- X[,accession.cols]
  X <- X[, -accession.cols]

  if (!use0) {
    X[X == 0] <- NA
  }
  
  res <- NULL
  for (i in 1:length(levels(group))) {
    
    X_ <- X[, group == levels(group)[i]]
    
    FUN <- switch(method, 
                  mean  = rowMeans, 
                  sum = rowSums, 
                  median = robustbase::rowMedians)
    
    X_ <- as.matrix(X_)
    
    res_tmp <- FUN(X_, na.rm = TRUE)
    
    if (!use0) {
      missingx <- apply(X_, 1, function(x) mean(is.na(x)))
      res_tmp[missingx > missing.limit | missingx == 1] <- NA
    }
    
    res <- cbind(res, res_tmp)
  }
  
  res <- as.data.frame(res)
  colnames(res) <- levels(group)
  res <- cbind(acc, res)
  return(res)
}

  

#### calculates peptide ratios for pairwise comparisons of groups
## D: dataset
## X: group1 (column name)
## Y: group2 (column name)
## useNA: if TRUE, results 0 and Inf are possible, otherwise ratio is NA if value for X or Y is NA
### result: fold changes (Y/X)
foldChange <- function(D, X, Y, useNA = FALSE) {
  FC <- D[, Y] / D[,X]
  
  if(useNA) {
    FC[is.na(D[,Y]) & !is.na(D[,X])] <- 0
    FC[is.na(D[,X]) & !is.na(D[,Y])] <- Inf
  }
  
  return(FC)
}







