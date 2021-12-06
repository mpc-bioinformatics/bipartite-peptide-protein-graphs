library(igraph)

source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")


### calculate isomorphism classes

################################################################################
#### D1_quant (without isoforms)

load("data/D1/D1_quant/preprocessed/Submatrices.RData")

for (i in 1:4) {
  for (j in (i+1):5) {
    print(paste0(i,j))
    S <- get(paste("Submatrix", i, j, sep = "_"))
    assign(paste("isomorph", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(isomorph_1_2,
     isomorph_1_3,
     isomorph_1_4,
     isomorph_1_5,
     isomorph_2_3,
     isomorph_2_4,
     isomorph_2_5,
     isomorph_3_4,
     isomorph_3_5,
     isomorph_4_5,
     file = "data/D1/D1_quant/isomorph_classes/isomorph_classes.RData")


################################################################################
##### for collaped peptide nodes:

load("data/D1/D1_quant/preprocessed/Submatrices_merged_Peptides.RData")

for (i in 1:4) {
  for (j in (i+1):5) {
    print(paste0(i,j))
    S <- get(paste("Submatrix_merged_Peptides", i, j, sep = "_"))
    assign(paste("isomorph_merged_Peptides", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(isomorph_merged_Peptides_1_2,
     isomorph_merged_Peptides_1_3,
     isomorph_merged_Peptides_1_4,
     isomorph_merged_Peptides_1_5,
     isomorph_merged_Peptides_2_3,
     isomorph_merged_Peptides_2_4,
     isomorph_merged_Peptides_2_5,
     isomorph_merged_Peptides_3_4,
     isomorph_merged_Peptides_3_5,
     isomorph_merged_Peptides_4_5,
     file = "data/D1/D1_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")


################################################################################
### complete isomorph list containing all 10 comparisons

load("data/D1/D1_quant/preprocessed/Submatrices.RData")
load("data/D1/D1_quant/preprocessed/Submatrices_merged_Peptides.RData")

Submatrix_all <- do.call(c, args = list(Submatrix_1_2,
                                        Submatrix_1_3,
                                        Submatrix_1_4,
                                        Submatrix_1_5,
                                        Submatrix_2_3,
                                        Submatrix_2_4,
                                        Submatrix_2_5,
                                        Submatrix_3_4,
                                        Submatrix_3_5,
                                        Submatrix_4_5
))


isomorph_all <- calculateIsomorphList(Submatrix_all)
save(isomorph_all, file = "data/D1/D1_quant/isomorph_classes/isomorph_classes_all.RData")


Submatrix_all_merged_Peptides <- do.call(c, args = list(Submatrix_merged_Peptides_1_2,
                                                        Submatrix_merged_Peptides_1_3,
                                                        Submatrix_merged_Peptides_1_4,
                                                        Submatrix_merged_Peptides_1_5,
                                                        Submatrix_merged_Peptides_2_3,
                                                        Submatrix_merged_Peptides_2_4,
                                                        Submatrix_merged_Peptides_2_5,
                                                        Submatrix_merged_Peptides_3_4,
                                                        Submatrix_merged_Peptides_3_5,
                                                        Submatrix_merged_Peptides_4_5
))

isomorph_all_merged_Peptides <- calculateIsomorphList(Submatrix_all_merged_Peptides)
save(isomorph_all_merged_Peptides, file = "data/D1/D1_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")

################################################################################
################################################################################
################################################################################
#### D2_quant (without isoforms)

load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices.RData")

for (i in 1:8) {
  for (j in (i+1):9) {
    S <- get(paste("Submatrix", i, j, sep = "_"))
    assign(paste("isomorph", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph"),
     file = "data/D2_without_isoforms/D2_quant/isomorph_classes/isomorph_classes.RData")

################################################################################
#### for collapsed peptide nodes:

load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices_merged_Peptides.RData")

for (i in 1:8) {
  for (j in (i+1):9) {
    S <- get(paste("Submatrix_merged_Peptides", i, j, sep = "_"))
    assign(paste("isomorph_merged_Peptides", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph_merged"),
     file = "data/D2_without_isoforms/D2_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")



################################################################################
#### complete isomorph list containing all comparisons

load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices.RData")
Submatrix_all <- do.call(c, args = mget(ls(pattern = "Submatrix")))

isomorph_all <- calculateIsomorphList(Submatrix_all)
save(isomorph_all, file = "data/D2_without_isoforms/D2_quant/isomorph_classes/isomorph_classes_all.RData")

load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices_merged_Peptides.RData")
Submatrix_all_merged_Peptides <- do.call(c, args = mget(ls(pattern = "Submatrix_merged_Peptides")))
isomorph_all_merged_Peptides <- calculateIsomorphList(Submatrix_all_merged_Peptides)
save(isomorph_all_merged_Peptides, file = "data/D2_without_isoforms/D2_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")



################################################################################
################################################################################
################################################################################
#### D3_quant (without isoforms)

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices.RData")

for (i in 1) {
  for (j in 2) {
    print(c(i,j))
    S <- get(paste("Submatrix", i, j, sep = "_"))
    assign(paste("isomorph", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph"),
     file = "data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes.RData")

################################################################################
#### for collapsed peptide nodes:

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")

for (i in 1) {
  for (j in 2) {
    S <- get(paste("Submatrix_merged_Peptides", i, j, sep = "_"))
    assign(paste("isomorph_merged_Peptides", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph_merged"),
     file = "data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")



################################################################################
#### complete isomorph list containing all comparisons

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices.RData")
Submatrix_all <- Submatrix_1_2

isomorph_all <- calculateIsomorphList(Submatrix_all)
save(isomorph_all, file = "data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes_all.RData")

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")
Submatrix_all_merged_Peptides <- Submatrix_merged_Peptides_1_2
isomorph_all_merged_Peptides <- calculateIsomorphList(Submatrix_all_merged_Peptides)
save(isomorph_all_merged_Peptides, file = "data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")




################################################################################
################################################################################
################################################################################
#### D3_quant (with isoforms)

load("data/D3/D3_quant/preprocessed/Submatrices.RData")

for (i in 1) {
  for (j in 2) {
    print(c(i,j))
    S <- get(paste("Submatrix", i, j, sep = "_"))
    assign(paste("isomorph", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph"),
     file = "data/D3/D3_quant/isomorph_classes/isomorph_classes.RData")

################################################################################
#### for collapsed peptide nodes:

load("data/D3/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")

for (i in 1) {
  for (j in 2) {
    S <- get(paste("Submatrix_merged_Peptides", i, j, sep = "_"))
    assign(paste("isomorph_merged_Peptides", i, j, sep = "_"), calculateIsomorphList(S))
  }
}

save(list =  ls(pattern = "isomorph_merged"),
     file = "data/D3/D3_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")



################################################################################
#### complete isomorph list containing all comparisons

load("data/D3/D3_quant/preprocessed/Submatrices.RData")
Submatrix_all <- do.call(c, args = mget(ls(pattern = "Submatrix")))

isomorph_all <- calculateIsomorphList(Submatrix_all)
save(isomorph_all, file = "data/D3/D3_quant/isomorph_classes/isomorph_classes_all.RData")

load("data/D3/D3_quant//preprocessed/Submatrices_merged_Peptides.RData")
Submatrix_all_merged_Peptides <- do.call(c, args = mget(ls(pattern = "Submatrix_merged_Peptides")))
isomorph_all_merged_Peptides <- calculateIsomorphList(Submatrix_all_merged_Peptides)
save(isomorph_all_merged_Peptides, file = "data/D3/D3_quant//isomorph_classes/isomorph_classes_all_merged_Peptides.RData")



