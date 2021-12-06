library(igraph)
library(openxlsx)

source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

################################################################################
#### calculate table with characteristics of all graphs ####

# Submatrix_merged_Peptides: Submatrix list with collapsed protein and peptide nodes
# Submatrices_full: Submatrix list without collapsed protein and peptide nodes
# fastalevel: TRUE -> fasta level, FALSE -> dataset level
# comparison: name of comparison, for quantitative level

calculate_subgraph_characteristics <- function(Submatrix_merged_Peptides, Submatrices_full,
                                               fastalevel = TRUE, comparison = NULL) {

  Data <- NULL

  for (i in 1:length(Submatrix_merged_Peptides)) {

    submatrix_full_tmp <- Submatrices_full[[i]]#[[1]]
    submatrix_tmp <- Submatrix_merged_Peptides[[i]]#[[1]]

    if ("list" %in% class(submatrix_tmp)) {          # class list, if fold changes are present (-> reach matrix with $X)
      submatrix_tmp <- submatrix_tmp$X
      submatrix_full_tmp <- submatrix_full_tmp$X
    } else   {
      submatrix_tmp <- submatrix_tmp
      submatrix_full_tmp <- submatrix_full_tmp
    }

    dim_submatrix_full_tmp <- dim(submatrix_full_tmp)
    dim_submatrix_tmp <- dim(submatrix_tmp)

    ### for the fastalevel, protein nodes are already collapsed, and must be seperated again
    if (!fastalevel) {
      dim_submatrix_full_tmp[2] <- length(strsplit(paste(colnames(submatrix_tmp), "", collapse = ";"), ";")[[1]])
    }

    subgraph_tmp <- igraph::convertToBipartiteGraph(submatrix_tmp)

    is_unique_peptide <- rowSums(submatrix_tmp) == 1
    nr_unique_pep_node_per_prot <- colSums(submatrix_tmp[is_unique_peptide,, drop = FALSE])
    nr_shared_pep_node_per_prot <- colSums(submatrix_tmp[!is_unique_peptide,, drop = FALSE])

    D_tmp <- data.frame(ID = i,
               Nr_prot_acc = dim_submatrix_full_tmp[2],
               Nr_prot_node = dim_submatrix_tmp[2],
               Nr_prot_node_only_unique_pep = sum(nr_unique_pep_node_per_prot > 0 & nr_shared_pep_node_per_prot == 0),
               Nr_prot_node_uniq_and_shared_pep = sum(nr_unique_pep_node_per_prot > 0 & nr_shared_pep_node_per_prot > 0),
               Nr_prot_node_only_shared_pep = sum(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0),
               Nr_prot_node_with_unique_pep = sum(nr_unique_pep_node_per_prot > 0),

               ratio_prot_nodes_with_and_without_unique_pep = sum(nr_unique_pep_node_per_prot > 0) / sum(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0),
               ratio_prot_nodes_without_and_with_unique_pep = sum(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0) / sum(nr_unique_pep_node_per_prot > 0),
               ratio_pep_nodes_prot_nodes = dim_submatrix_tmp[1] / dim_submatrix_tmp[2],
               mean_pep_node_per_prot_node = mean(colSums(submatrix_tmp)),
               mean_unique_pep_node_per_prot_node = mean(nr_unique_pep_node_per_prot),
               mean_shared_pep_node_per_prot_node = mean(nr_shared_pep_node_per_prot),

               has_multiple_prot = (dim_submatrix_tmp[2] > 1),
               has_prot_without_unique_pep = any(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0),

               Nr_pep_seq = dim_submatrix_full_tmp[1],
               Nr_pep_node = dim_submatrix_tmp[1],
               Nr_pep_node_unique = sum(is_unique_peptide),
               Nr_pep_node_shared = sum(!is_unique_peptide),

               Nr_edge = gsize(subgraph_tmp)
    )

    if (!fastalevel) {
      D_tmp <- data.frame(comparison = comparison, D_tmp)
    }

    Data <- rbind(Data, D_tmp)
  }

  return(Data)
}


################################################################################
################################################################################
################################################################################
#### D1_fasta (without isoforms)

Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min5AA_fast.rds")
Submatrices_full <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min5AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min5AA.xlsx", keepNA = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min6AA_fast.rds")
Submatrices_full <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min6AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min6AA.xlsx", keepNA = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min7AA_fast.rds")
Submatrices_full <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min7AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min7AA.xlsx", keepNa = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_merged_Peptides_D1_fasta_min9AA_fast.rds")
Submatrices_full <- readRDS("data/D1/D1_fasta/preprocessed/Submatrix_D1_fasta_min9AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min9AA.xlsx", keepNa = TRUE)

################################################################################
######## D1_quant (without isoforms)

load("data/D1/D1_quant/preprocessed/Submatrices.RData")
load("data/D1/D1_quant/preprocessed/Submatrices_merged_Peptides.RData")
D2 <- NULL
for (i in 1:4) {
  for (j in (i+1):5) {
    comparison <- paste(i, j, sep = "_")

    S_merged_peptides <- get(paste0("Submatrix_merged_Peptides_", i, "_", j))
    S_full <- get(paste0("Submatrix_", i, "_", j))
    D_tmp2 <- calculate_subgraph_characteristics(Submatrix_merged_Peptides = S_merged_peptides,
                                                Submatrices_full = S_full,
                                                fastalevel = FALSE, comparison = comparison)
    D2 <- rbind(D2, D_tmp2)
  }
}
write.xlsx(D2, "data/D1/D1_quant/table_subgraph_characteristics_D1_quant.xlsx", keepNA = TRUE)



################################################################################
################################################################################
################################################################################
#### D2_fasta (without isoforms)

Submatrix_merged_Peptides <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min5AA_fast.rds")
Submatrices_full <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_D2_fasta_min5AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min5AA.xlsx", keepNA = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min6AA_fast.rds")
Submatrices_full <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_D2_fasta_min6AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min6AA.xlsx", keepNA = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min7AA_fast.rds")
Submatrices_full <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_D2_fasta_min7AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min7AA.xlsx", keepNa = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_merged_Peptides_D2_fasta_min9AA_fast.rds")
Submatrices_full <- readRDS("data/D2_without_isoforms/D2_fasta/preprocessed/Submatrix_D2_fasta_min9AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min9AA.xlsx", keepNa = TRUE)


################################################################################
######## D2_quant (without isoforms)

load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices.RData")
load("data/D2_without_isoforms/D2_quant/preprocessed/Submatrices_merged_Peptides.RData")
D2 <- NULL
for (i in 1:8) {
  for (j in (i+1):9) {
    comparison <- paste(i, j, sep = "_")

    S_merged_peptides <- get(paste0("Submatrix_merged_Peptides_", i, "_", j))
    S_full <- get(paste0("Submatrix_", i, "_", j))
    D_tmp2 <- calculate_subgraph_characteristics(Submatrix_merged_Peptides = S_merged_peptides,
                                                 Submatrices_full = S_full,
                                                 fastalevel = FALSE, comparison = comparison)
    D2 <- rbind(D2, D_tmp2)
  }
}
write.xlsx(D2, "data/D2_without_isoforms/D2_quant/table_subgraph_characteristics_D2_quant.xlsx", keepNA = TRUE)


################################################################################
################################################################################
################################################################################
#### D3_fasta (without isoforms)

Submatrix_merged_Peptides <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_without_isoforms_fasta_min5AA_fast.rds")
Submatrices_full <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_D3_without_isoforms_fasta_min5AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min5AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_without_isoforms_fasta_min6AA_fast.rds")
Submatrices_full <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_D3_without_isoforms_fasta_min6AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min6AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_without_isoforms_fasta_min7AA_fast.rds")
Submatrices_full <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_D3_without_isoforms_fasta_min7AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min7AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_without_isoforms_fasta_min9AA_fast.rds")
Submatrices_full <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/Submatrix_D3_without_isoforms_fasta_min9AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min9AA.xlsx", keepNA = TRUE, overwrite = TRUE)



################################################################################
######## D3_quant (without isoforms)

load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices.RData")
load("data/D3_without_isoforms/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")
D3 <- NULL
for (i in 1) {
  for (j in 2) {
    comparison <- paste(i, j, sep = "_")

    S_merged_peptides <- get(paste0("Submatrix_merged_Peptides_", i, "_", j))
    S_full <- get(paste0("Submatrix_", i, "_", j))
    D_tmp2 <- calculate_subgraph_characteristics(Submatrix_merged_Peptides = S_merged_peptides,
                                                 Submatrices_full = S_full,
                                                 fastalevel = FALSE, comparison = comparison)
    D3 <- rbind(D3, D_tmp2)
  }
}
write.xlsx(D3, "data/D3_without_isoforms/D3_quant/table_subgraph_characteristics_D3_quant.xlsx", keepNA = TRUE, overwrite = TRUE)


################################################################################
################################################################################
################################################################################
#### D3_fasta (with isoforms)

Submatrix_merged_Peptides <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_fasta_min5AA_fast.rds")
Submatrices_full <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_D3_fasta_min5AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min5AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_fasta_min6AA_fast.rds")
Submatrices_full <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_D3_fasta_min6AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min6AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_fasta_min7AA_fast.rds")
Submatrices_full <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_D3_fasta_min7AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min7AA.xlsx", keepNA = TRUE, overwrite = TRUE)

Submatrix_merged_Peptides <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_merged_Peptides_D3_fasta_min9AA_fast.rds")
Submatrices_full <- readRDS("data/D3/D3_fasta/preprocessed/Submatrix_D3_fasta_min9AA.rds")
D <- calculate_subgraph_characteristics(Submatrix_merged_Peptides, Submatrices_full)
write.xlsx(D, "data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min9AA.xlsx", keepNA = TRUE, overwrite = TRUE)



################################################################################
######## D3_quant (with isoforms)

load("data/D3/D3_quant/preprocessed/Submatrices.RData")
load("data/D3/D3_quant/preprocessed/Submatrices_merged_Peptides.RData")
D3 <- NULL
for (i in 1) {
  for (j in 2) {
    comparison <- paste(i, j, sep = "_")

    S_merged_peptides <- get(paste0("Submatrix_merged_Peptides_", i, "_", j))
    S_full <- get(paste0("Submatrix_", i, "_", j))
    D_tmp2 <- calculate_subgraph_characteristics(Submatrix_merged_Peptides = S_merged_peptides,
                                                 Submatrices_full = S_full,
                                                 fastalevel = FALSE, comparison = comparison)
    D3 <- rbind(D3, D_tmp2)
  }
}
write.xlsx(D3, "data/D3/D3_quant/table_subgraph_characteristics_D3_quant.xlsx", keepNA = TRUE, overwrite = TRUE)

