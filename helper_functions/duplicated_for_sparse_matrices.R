#### duplicated function for class dgCMatrix

duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
  MARGIN <- as.integer(MARGIN)
  n <- nrow(dgCMat)
  p <- ncol(dgCMat)
  J <- rep(1:p, diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    RowLst <- split(x, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    ColLst <- split(x, J)
    is_empty <- setdiff(1:p, J)
    result <- duplicated.default(ColLst)
  } else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  
  if(any(is_empty)){
    out <- logical(if(MARGIN == 1L) n else p)
    out[-is_empty] <- result
    if(length(is_empty) > 1)
      out[is_empty[-1]] <- TRUE
    result <- out
  }
  
  result
}



# taken from this stackoverflow post by JZL003:
# https://stackoverflow.com/questions/51457184/is-there-a-method-with-r-function-duplicated-for-a-sparse-matrix-from-matrix
# under the CC BY-SA 4.0 license (https://creativecommons.org/licenses/by-sa/4.0/)
