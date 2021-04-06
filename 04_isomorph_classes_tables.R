library(openxlsx)


################################################################################
#### calculate tables with information of on isomorph classes ####

###############################################################################
#### get different characteristics/attributes from a graph
## G: graph (igraoh object)
getInfoFromGraph <- function(G) {
  require(igraph)
  
  nr_node <- length(V(G))
  nr_node_protein <- sum(V(G)$type)
  nr_node_peptide <- sum(!V(G)$type)
  
  degrees_peptide_node <- unname(degree(G))[!V(G)$type]
  nr_node_unique_peptides <- sum(degrees_peptide_node == 1)
  nr_node_shared_peptides <- nr_node_peptide - nr_node_unique_peptides
  
  nr_edges <- gsize(G)
  
  return(c(nr_node, nr_node_protein, nr_node_peptide, nr_node_unique_peptides, 
           nr_node_shared_peptides, nr_edges))
  
}


################################################################################
#### function to generate table:

tab_isomorph_classes <- function (isomorph, file) {
  tab <- data.frame(id = 1:(length(isomorph$isomorph_list)),
                    nr_nodes = rep(NA, length(isomorph$isomorph_list)), 
                    nr_nodes_protein = rep(NA, length(isomorph$isomorph_list)), 
                    nr_nodes_peptide = rep(NA, length(isomorph$isomorph_list)), 
                    nr_nodes_unique_peptide = rep(NA, length(isomorph$isomorph_list)), 
                    nr_nodes_shared_peptide = rep(NA, length(isomorph$isomorph_list)), 
                    nr_edges = rep(NA, length(isomorph$isomorph_list)), 
                    nr_occurrences = rep(NA, length(isomorph$isomorph_list))) 
  
  for (i in 1:length(isomorph$isomorph_list)) {
    ind <- isomorph$isomorph_list[[i]][1]
    G <- isomorph$Graphs[[ind]]
    nr_occurrences <- length(isomorph$isomorph_list[[i]])
    res <- c(i, getInfoFromGraph(G), nr_occurrences)
    tab[i,] <- res
  }
  write.xlsx(tab, file)
}

################################################################################
#### D1_fasta_min5AA ####
load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min5AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D1/D1_fasta/table_isomorph_classes_D1_fasta_min5AA.xlsx")

#### D1_fasta_min6AA ####
load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min6AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D1/D1_fasta/table_isomorph_classes_D1_fasta_min6AA.xlsx")

#### D1_fasta_min7AA ####
load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D1/D1_fasta/table_isomorph_classes_D1_fasta_min7AA.xlsx")

#### D1_fasta_min9AA ####
load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min9AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D1/D1_fasta/table_isomorph_classes_D1_fasta_min9AA.xlsx")

#### D1_quant ####
load("data/D1/D1_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
tab_isomorph_classes(isomorph_all_merged_Peptides, file = "data/D1/D1_quant/table_isomorph_classes_D1_quant.xlsx")


################################################################################
#### D2_fasta_min5AA ####
load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min5AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D2/D2_fasta/table_isomorph_classes_D2_fasta_min5AA.xlsx")

#### D2_fasta_min6AA ####
load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D2/D2_fasta/table_isomorph_classes_D2_fasta_min6AA.xlsx")

#### D2_fasta_min7AA ####
load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min7AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D2/D2_fasta/table_isomorph_classes_D2_fasta_min7AA.xlsx")

#### D2_fasta_min9AA ####
load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min9AA_fast.RData")
tab_isomorph_classes(isomorph, file = "data/D2/D2_fasta/table_isomorph_classes_D2_fasta_min9AA.xlsx")

#### D2_quant ####
load("data/D2/D2_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
tab_isomorph_classes(isomorph_all_merged_Peptides, file = "data/D2/D2_quant/table_isomorph_classes_D2_quant.xlsx")


