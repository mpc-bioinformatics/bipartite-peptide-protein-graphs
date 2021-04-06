library(igraph)    # for functionality regarding (bipartite graphs)
library(pbapply)   # for progress bars for apply functions

source("helper_functions/calculate_biadjacency_matrix_and_submatrices.R")
source("helper_functions/duplicated_for_sparse_matrices.R")
source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")


################################################################################
#### calculate connected components (subgraphs) from biadjacency matrix ####

################################################################################
#### D1_fasta_min5AA

matrix_01 <- readRDS("data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.rds")

ind <- which(colSums(matrix_01) == 0)  # proteins that do not have peptides anymore
matrix_01 <- matrix_01[,-ind]
saveRDS(matrix_01, "data/D1/D1_fasta/preprocessed/01_matrix_D1_fasta_min5AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min5AA.rds")
saveRDS(subgraphs, file = "data/D1/D1_fasta/preprocessed/Subgraphs_D1_fasta_min5AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min5AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = TRUE)
saveRDS(S_proteingroups, file = "data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min5AA_fast.rds")

#### collapse peptide nodes
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min5AA_fast.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min5AA_fast.rds")

#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min5AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min5AA_fast.RData")




################################################################################
#### D1_fasta_min6AA

matrix_01 <- readRDS("data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.rds")

ind <- which(nchar(rownames(matrix_01)) < 6)  # too small peptides
matrix_01 <- matrix_01[-ind,]
ind2 <- which(colSums(matrix_01) == 0)  # proteins that do not have peptides anymore
matrix_01 <- matrix_01[,-ind2]
saveRDS(matrix_01, "data/D1/D1_fasta/preprocessed/01_matrix_D1_fasta_min6AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min6AA.rds")
saveRDS(subgraphs, file = "data/D1/D1_fasta/preprocessed/Subgraphs_D1_fasta_min6AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min6AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = TRUE)
saveRDS(S_proteingroups, file = "data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min6AA_fast.rds")

#### collapse peptide nodes
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min6AA_fast.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min6AA_fast.rds")


#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min6AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min6AA_fast.RData")


################################################################################
#### D1_fasta_min7AA

matrix_01 <- readRDS("data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.rds")

ind <- which(nchar(rownames(matrix_01)) < 7)
matrix_01 <- matrix_01[-ind,]
ind2 <- which(colSums(matrix_01) == 0)  # proteins that do not have peptides anymore
matrix_01 <- matrix_01[,-ind2]
saveRDS(matrix_01, "data/D1/D1_fasta/preprocessed/01_matrix_D1_fasta_min7AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min7AA.rds")
saveRDS(subgraphs, file = "data/D1/D1_fasta/preprocessed/Subgraphs_D1_fasta_min7AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min7AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = FALSE)
saveRDS(S_proteingroups, file = "data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min7AA.rds")

#### collapse peptide nodes
S <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_proteingroups_D1_fasta_min7AA.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min7AA_fast.rds")

#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min7AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData")




################################################################################
################################################################################
################################################################################
#### D2_fasta_min5AA

matrix_01 <- readRDS("data/D2/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.rds")

ind <- which(colSums(matrix_01) == 0)
matrix_01 <- matrix_01[,-ind]
saveRDS(matrix_01, "data/D2/D2_fasta/preprocessed/01_matrix_D2_fasta_min5AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min5AA.rds")
saveRDS(subgraphs, file = "data/D2/D2_fasta/preprocessed/Subgraphs_D2_fasta_min5AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min5AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = TRUE)
saveRDS(S_proteingroups, file = "data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min5AA_fast.rds")

#### collapse peptide nodes
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min5AA_fast.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min5AA_fast.rds")

#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min5AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min5AA_fast.RData")



################################################################################
#### D2_fasta_min6AA

matrix_01 <- readRDS("data/D2/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.rds")

ind <- which(nchar(rownames(matrix_01)) < 6)
matrix_01 <- matrix_01[-ind,]
ind2 <- which(colSums(matrix_01) == 0)
matrix_01 <- matrix_01[,-ind2]
saveRDS(matrix_01, "data/D2/D2_fasta/preprocessed/01_matrix_D2_fasta_min6AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min6AA.rds")
saveRDS(subgraphs, file = "data/D2/D2_fasta/preprocessed/Subgraphs_D2_fasta_min6AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min6AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = FALSE)
saveRDS(S_proteingroups, file = "data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min6AA.rds")

#### collapse peptide nodes
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min6AA.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min6AA_fast.rds")

#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min6AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData")


################################################################################
#### D2_fasta_min7AA

matrix_01 <- readRDS("data/D2/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.rds")

ind <- which(nchar(rownames(matrix_01)) < 7)
matrix_01 <- matrix_01[-ind,]
ind2 <- which(colSums(matrix_01) == 0)
matrix_01 <- matrix_01[,-ind2]
saveRDS(matrix_01, "data/D2/D2_fasta/preprocessed/01_matrix_D2_fasta_min7AA.rds")

G <- igraph::graph_from_incidence_matrix(matrix_01) # calculate complete bipartite graph
subgraphs <- igraph::decompose(G)                   # calculate connected components        

### back to sub-adjacency matrices
Submatrix <- pblapply(subgraphs, as_incidence_matrix, sparse = TRUE)
Submatrix <- pblapply(Submatrix, as, "dgCMatrix")

saveRDS(Submatrix, file = "data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min7AA.rds")
saveRDS(subgraphs, file = "data/D2/D2_fasta/preprocessed/Subgraphs_D2_fasta_min7AA.rds")


#### collapse protein nodes:
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_D2_fasta_min7AA.rds")
S_proteingroups <- form_proteingroups(S, sparse = TRUE, fast = TRUE)
saveRDS(S_proteingroups, file = "data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min7AA_fast.rds")

#### collapse peptide nodes
S <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_proteingroups_D2_fasta_min7AA_fast.rds")
Submatrix_merged_Peptides <- merge_Peptides(S, sparse = TRUE, fc = FALSE, fast = TRUE)
saveRDS(Submatrix_merged_Peptides, file = "data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min7AA_fast.rds")

#### calculate isomorphism classes
Submatrix_merged_Peptides <- readRDS("data/D2/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min7AA_fast.rds")
isomorph <- calculateIsomorphList(Submatrix_merged_Peptides) 
save(isomorph, file = "data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min7AA_fast.RData")








