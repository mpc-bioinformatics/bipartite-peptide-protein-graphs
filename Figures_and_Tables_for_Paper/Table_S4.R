library(openxlsx)
library(xtable)

D3_fasta <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min7AA.xlsx")
isomorphlist_D3_fasta <- get(load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min7AA_fast.RData"))
D3_quant <- read.xlsx("data/D3_without_isoforms/D3_quant/table_subgraph_characteristics_D3_quant.xlsx")
isomorphlist_D3_quant <- get(load("data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData"))

D3_iso_fasta <- read.xlsx("data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min7AA.xlsx")
isomorphlist_D3_iso_fasta <- get(load("data/D3/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_fasta_min7AA_fast.RData"))
D3_iso_quant <- read.xlsx("data/D3/D3_quant/table_subgraph_characteristics_D3_quant.xlsx")
isomorphlist_D3_iso_quant <- get(load("data/D3/D3_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData"))


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


################################################################################
### Table S4: Comparison of D3 and D3_iso

RES <- mapply(table1_calculate_numbers, D = list(D3_fasta, D3_iso_fasta, D3_quant, D3_iso_quant),
              isomorph = list(isomorphlist_D3_fasta, isomorphlist_D3_iso_fasta, isomorphlist_D3_quant, isomorphlist_D3_iso_quant))
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")

RES2 <- RES
colnames(RES2) <- c("D3_fasta", "D3_iso_fasta", "D3_quant", "D3_iso_quant")
write.xlsx(RES2, "Paper/Paper 1/tables/TableS1_Paper.xlsx", keepNA = TRUE,
           row.names = TRUE, overwrite = TRUE)

xtable(RES, digits = 0)
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Mon Dec 20 13:54:05 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrr}
# \hline
# & 1 & 2 & 3 & 4 \\
# \hline
# proteins & 81548 & 103541 & 17585 & 22216 \\
# protein.nodes & 80856 & 102463 & 10969 & 11672 \\
# peptides & 3050340 & 3154386 & 30369 & 29490 \\
# peptide.nodes & 148555 & 184932 & 10540 & 10895 \\
# edges & 431391 & 647037 & 23802 & 27038 \\
# graphs & 20270 & 20048 & 5267 & 5162 \\
# Graphs.with.1.protein.node & 10129 & 8948 & 3315 & 3106 \\
# isomorphism.classes & 2305 & 5657 & 459 & 537 \\
# protein.nodes.1 & 6472 & 8895 & 57 & 58 \\
# peptide.nodes.1 & 14993 & 19123 & 65 & 67 \\
# edges.1 & 56950 & 87852 & 350 & 359 \\
# protein.nodes.2 & 306 & 434 & 35 & 37 \\
# peptide.nodes.2 & 757 & 940 & 34 & 34 \\
# edges.2 & 2884 & 4613 & 181 & 185 \\
# \hline
# \end{tabular}
# \end{table}

