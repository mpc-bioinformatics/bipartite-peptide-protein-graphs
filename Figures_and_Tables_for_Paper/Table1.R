library(openxlsx)
library(xtable)

################################################################################
### Table 1: impact of peptide length

################################################################################
### read in tables and data:
D1_fasta_5 <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min5AA.xlsx")
D1_fasta_6 <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min6AA.xlsx")
D1_fasta_7 <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min7AA.xlsx")
D1_fasta_9 <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min9AA.xlsx")

isomorphlist_D1_5 <- get(load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min5AA_fast.RData"))
isomorphlist_D1_6 <- get(load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min6AA_fast.RData"))
isomorphlist_D1_7 <- get(load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData"))
isomorphlist_D1_9 <- get(load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min9AA_fast.RData"))

D2_fasta_5 <- read.xlsx("data/D2/D2_fasta/table_subgraph_characteristics_D2_fasta_min5AA.xlsx")
D2_fasta_6 <- read.xlsx("data/D2/D2_fasta/table_subgraph_characteristics_D2_fasta_min6AA.xlsx")
D2_fasta_7 <- read.xlsx("data/D2/D2_fasta/table_subgraph_characteristics_D2_fasta_min7AA.xlsx")
D2_fasta_9 <- read.xlsx("data/D2/D2_fasta/table_subgraph_characteristics_D2_fasta_min9AA.xlsx")

isomorphlist_D2_5 <- get(load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min5AA_fast.RData"))
isomorphlist_D2_6 <- get(load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData"))
isomorphlist_D2_7 <- get(load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min7AA_fast.RData"))
isomorphlist_D2_9 <- get(load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min9AA_fast.RData"))


################################################################################
### function to calculate table:

table1_calculate_numbers <- function(D, isomorph) {
  
  ind_zero_pep <- which(D$Nr_pep_node == 0)
  if (length(ind_zero_pep) > 0) D <- D[-ind_zero_pep, ]
  
  ### largest graph (in terms of number of protein nodes)
  ind_largest <- which.max(D$Nr_prot_node)
  
  ### 2nd largest graph (in terms of number of protein nodes)
  ind_largest2 <- which(D$Nr_prot_node == sort(D$Nr_prot_node, decreasing = TRUE)[2])
  
  c( 
    sum(D$Nr_prot_acc),                 ### Nr of proteins
    sum(D$Nr_prot_node),                ### Nr of protein nodes
    sum(D$Nr_pep_seq),                  ### Nr of peptides
    sum(D$Nr_pep_node),                 ### Nr of peptide nodes
    sum(D$Nr_edge),                     ### Nr of edges
    
    as.integer(nrow(D)),                ### Nr of graphs
    as.integer(sum(D$Nr_prot_node == 1 & D$Nr_pep_node == 1)),                 ### Nr of graphs with only 1 protein node
    as.integer(length(isomorph$isomorph_list) - 1*(length(ind_zero_pep) > 0)), ### Nr of isomorphism classes

    D$Nr_prot_node[ind_largest],        ### largest system
    D$Nr_pep_node[ind_largest],
    D$Nr_edge[ind_largest],
    
    D$Nr_prot_node[ind_largest2],       ### 2nd largest system
    D$Nr_pep_node[ind_largest2],
    D$Nr_edge[ind_largest2]
  
  )
  
}


#######################################################################################
### compile table:


RES <- mapply(table1_calculate_numbers, D = list(D1_fasta_5, D1_fasta_6, D1_fasta_7, D1_fasta_9, 
                                                 D2_fasta_5, D2_fasta_6, D2_fasta_7, D2_fasta_9), 
       isomorph = list(isomorphlist_D1_5, isomorphlist_D1_6, isomorphlist_D1_7, isomorphlist_D1_9,
                       isomorphlist_D2_5, isomorphlist_D2_6, isomorphlist_D2_7, isomorphlist_D2_9)) 
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs", 
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges", 
                   "protein nodes", "peptide nodes", "edges")

xtable(RES, digits = 0)

