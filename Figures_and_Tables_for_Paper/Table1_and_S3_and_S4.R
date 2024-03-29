library(openxlsx)
library(xtable)

################################################################################
### Table 1, S3 and S4: impact of peptide length for the database level

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

D2_fasta_5 <- read.xlsx("data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min5AA.xlsx")
D2_fasta_6 <- read.xlsx("data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min6AA.xlsx")
D2_fasta_7 <- read.xlsx("data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min7AA.xlsx")
D2_fasta_9 <- read.xlsx("data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min9AA.xlsx")
isomorphlist_D2_5 <- get(load("data/D2_without_isoforms/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min5AA_fast.RData"))
isomorphlist_D2_6 <- get(load("data/D2_without_isoforms/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData"))
isomorphlist_D2_7 <- get(load("data/D2_without_isoforms/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min7AA_fast.RData"))
isomorphlist_D2_9 <- get(load("data/D2_without_isoforms/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min9AA_fast.RData"))

D3_fasta_5 <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min5AA.xlsx")
D3_fasta_6 <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min6AA.xlsx")
D3_fasta_7 <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min7AA.xlsx")
D3_fasta_9 <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min9AA.xlsx")
isomorphlist_D3_5 <- get(load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min5AA_fast.RData"))
isomorphlist_D3_6 <- get(load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min6AA_fast.RData"))
isomorphlist_D3_7 <- get(load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min7AA_fast.RData"))
isomorphlist_D3_9 <- get(load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min9AA_fast.RData"))


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
### Table 1: D1, min 5, 6, 7 and 9 AA

RES <- mapply(table1_calculate_numbers, D = list(D1_fasta_5, D1_fasta_6, D1_fasta_7, D1_fasta_9),
       isomorph = list(isomorphlist_D1_5, isomorphlist_D1_6, isomorphlist_D1_7, isomorphlist_D1_9))
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")

RES2 <- RES
colnames(RES2) <- c("D1_fasta_5", "D1_fasta_6", "D1_fasta_7", "D1_fasta_9")
write.xlsx(RES2, "Paper/Paper 1/tables/Table1.xlsx", keepNA = TRUE,
           row.names = TRUE, overwrite = TRUE)

xtable(RES, digits = 0)
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Thu Dec 02 14:13:13 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrr}
# \hline
# & 1 & 2 & 3 & 4 \\
# \hline
# proteins & 52760 & 52752 & 52741 & 52679 \\
# protein.nodes & 52405 & 52388 & 52364 & 52250 \\
# peptides & 2872496 & 2777410 & 2642464 & 2365974 \\
# peptide.nodes & 133137 & 103517 & 96241 & 92346 \\
# edges & 444272 & 277656 & 241601 & 221186 \\
# graphs & 2916 & 10356 & 15402 & 17114 \\
# Graphs.with.1.protein.node & 2004 & 4809 & 6355 & 7147 \\
# isomorphism.classes & 206 & 1474 & 2583 & 2668 \\
# protein.nodes.1 & 47525 & 25006 & 1632 & 939 \\
# peptide.nodes.1 & 126315 & 57763 & 3971 & 2415 \\
# edges.1 & 433174 & 185147 & 14093 & 8751 \\
# protein.nodes.2 & 26 & 91 & 1104 & 529 \\
# peptide.nodes.2 & 75 & 177 & 2940 & 1046 \\
# edges.2 & 288 & 486 & 10541 & 2484 \\
# \hline
# \end{tabular}
# \end{table}


#######################################################################################
### Table S3: D2, min 5, 6, 7 and 9 AA

RES <- mapply(table1_calculate_numbers, D = list(D2_fasta_5, D2_fasta_6, D2_fasta_7, D2_fasta_9),
              isomorph = list(isomorphlist_D2_5, isomorphlist_D2_6, isomorphlist_D2_7, isomorphlist_D2_9))
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")

RES2 <- RES
colnames(RES2) <- c("D2_fasta_5", "D2_fasta_6", "D2_fasta_7", "D2_fasta_9")
write.xlsx(RES2, "Paper/Paper 1/tables/Table2.xlsx", keepNA = TRUE, row.names = TRUE, overwrite = TRUE)

xtable(RES, digits = 0)
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Thu Dec 02 14:18:25 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrr}
# \hline
# & 1 & 2 & 3 & 4 \\
# \hline
# proteins & 6336 & 6335 & 6333 & 6333 \\
# protein.nodes & 6265 & 6264 & 6263 & 6263 \\
# peptides & 752532 & 718476 & 679995 & 603418 \\
# peptide.nodes & 12962 & 8132 & 7563 & 7366 \\
# edges & 26792 & 15029 & 13653 & 12677 \\
# graphs & 1649 & 4908 & 5471 & 5604 \\
# Graphs.with.1.protein.node & 1506 & 4296 & 5045 & 5264 \\
# isomorphism.classes & 18 & 64 & 41 & 36 \\
# protein.nodes.1 & 4413 & 116 & 69 & 69 \\
# peptide.nodes.1 & 10907 & 356 & 264 & 253 \\
# edges.1 & 24434 & 1292 & 3329 & 3033 \\
# protein.nodes.2 & 21 & 70 & 52 & 50 \\
# peptide.nodes.2 & 37 & 267 & 203 & 181 \\
# edges.2 & 131 & 3359 & 743 & 613 \\
# \hline
# \end{tabular}
# \end{table}


####################################################################################
#### Table S4: D3, min 5, 6, 7 and 9 AA

RES <- mapply(table1_calculate_numbers, D = list(D3_fasta_5, D3_fasta_6, D3_fasta_7, D3_fasta_9),
              isomorph = list(isomorphlist_D3_5, isomorphlist_D3_6, isomorphlist_D3_7, isomorphlist_D3_9))
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")

RES2 <- RES
colnames(RES2) <- c("D3_fasta_5", "D3_fasta_6", "D3_fasta_7", "D3_fasta_9")
write.xlsx(as.data.frame(RES2), "Paper/Paper 1/tables/Table3.xlsx",
           row.names = TRUE, overwrite = TRUE)

xtable(RES, digits = 0)
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Mon Dec 20 14:28:33 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrr}
# \hline
# & 1 & 2 & 3 & 4 \\
# \hline
# proteins & 81591 & 81572 & 81548 & 81440 \\
# protein.nodes & 80932 & 80897 & 80856 & 80676 \\
# peptides & 3309331 & 3204104 & 3050340 & 2733226 \\
# peptide.nodes & 192162 & 157556 & 148555 & 143264 \\
# edges & 699493 & 480830 & 431391 & 401884 \\
# graphs & 4576 & 14177 & 20270 & 22327 \\
# Graphs.with.1.protein.node & 3722 & 8178 & 10129 & 11088 \\
# isomorphism.classes & 253 & 2305 & 4198 & 4522 \\
# protein.nodes.1 & 74157 & 40266 & 6472 & 2203 \\
# peptide.nodes.1 & 183383 & 89697 & 14993 & 5106 \\
# edges.1 & 685438 & 315173 & 56950 & 17897 \\
# protein.nodes.2 & 27 & 86 & 306 & 229 \\
# peptide.nodes.2 & 51 & 126 & 757 & 454 \\
# edges.2 & 198 & 1412 & 2884 & 2084 \\
# \hline
# \end{tabular}
# \end{table}
