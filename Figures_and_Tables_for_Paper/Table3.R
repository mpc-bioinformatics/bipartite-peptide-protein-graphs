library(openxlsx)
library(xtable)
library(tidyverse)
library(matrixStats)



################################################################################
### function to calculate table:

table3_calculate_numbers <- function(D, isomorph) {

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
################################################################################
################################################################################
### dataset D1

D1_tab <- read.xlsx("data/D1/D1_quant/table_subgraph_characteristics_D1_quant.xlsx")

load("data/D1/D1_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")


RES <- mapply(table3_calculate_numbers, D = list(filter(D1_tab, comparison == "1_2"),
                                                 filter(D1_tab, comparison == "1_3"),
                                                 filter(D1_tab, comparison == "1_4"),
                                                 filter(D1_tab, comparison == "1_5"),
                                                 filter(D1_tab, comparison == "2_3"),
                                                 filter(D1_tab, comparison == "2_4"),
                                                 filter(D1_tab, comparison == "2_5"),
                                                 filter(D1_tab, comparison == "3_4"),
                                                 filter(D1_tab, comparison == "3_5"),
                                                 filter(D1_tab, comparison == "4_5")),
              isomorph = list(isomorph_merged_Peptides_1_2,
                              isomorph_merged_Peptides_1_3,
                              isomorph_merged_Peptides_1_4,
                              isomorph_merged_Peptides_1_5,
                              isomorph_merged_Peptides_2_3,
                              isomorph_merged_Peptides_2_4,
                              isomorph_merged_Peptides_2_5,
                              isomorph_merged_Peptides_3_4,
                              isomorph_merged_Peptides_3_5,
                              isomorph_merged_Peptides_4_5))
rownames(RES) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")

RES2 <- RES
colnames(RES2) <- c("1_2", "1_3", "1_4", "1_5", "2_3", "2_4", "2_5", "3_4", "3_5", "4_5")
write.xlsx(RES2, "data/D1/D1_quant/draft_Table3_Paper.xlsx", keepNA = TRUE, row.names = TRUE)


RES3 <- cbind(round(rowMeans(RES), 0), matrixStats::rowMins(RES), matrixStats::rowMaxs(RES))
colnames(RES3) <- c("mean", "min", "max")


################################################################################
################################################################################
################################################################################
#### dataset D2


D2_tab <- read.xlsx("data/D2/D2_quant/table_subgraph_characteristics_D2_quant.xlsx")

load("data/D2/D2_quant/isomorph_classes/isomorph_classes_merged_Peptides.RData")


RES_2 <- mapply(table3_calculate_numbers, D = list(filter(D2_tab, comparison == "1_2"),
                                                   filter(D2_tab, comparison == "1_3"),
                                                   filter(D2_tab, comparison == "1_4"),
                                                   filter(D2_tab, comparison == "1_5"),
                                                   filter(D2_tab, comparison == "1_6"),
                                                   filter(D2_tab, comparison == "1_7"),
                                                   filter(D2_tab, comparison == "1_8"),
                                                   filter(D2_tab, comparison == "1_9"),
                                                   filter(D2_tab, comparison == "2_3"),
                                                   filter(D2_tab, comparison == "2_4"),
                                                   filter(D2_tab, comparison == "2_5"),
                                                   filter(D2_tab, comparison == "2_6"),
                                                   filter(D2_tab, comparison == "2_7"),
                                                   filter(D2_tab, comparison == "2_8"),
                                                   filter(D2_tab, comparison == "2_9"),
                                                   filter(D2_tab, comparison == "3_4"),
                                                   filter(D2_tab, comparison == "3_5"),
                                                   filter(D2_tab, comparison == "3_6"),
                                                   filter(D2_tab, comparison == "3_7"),
                                                   filter(D2_tab, comparison == "3_8"),
                                                   filter(D2_tab, comparison == "3_9"),
                                                   filter(D2_tab, comparison == "4_5"),
                                                   filter(D2_tab, comparison == "4_6"),
                                                   filter(D2_tab, comparison == "4_7"),
                                                   filter(D2_tab, comparison == "4_8"),
                                                   filter(D2_tab, comparison == "4_9"),
                                                   filter(D2_tab, comparison == "5_6"),
                                                   filter(D2_tab, comparison == "5_7"),
                                                   filter(D2_tab, comparison == "5_8"),
                                                   filter(D2_tab, comparison == "5_9"),
                                                   filter(D2_tab, comparison == "6_7"),
                                                   filter(D2_tab, comparison == "6_8"),
                                                   filter(D2_tab, comparison == "6_9"),
                                                   filter(D2_tab, comparison == "7_8"),
                                                   filter(D2_tab, comparison == "7_9"),
                                                   filter(D2_tab, comparison == "8_9")),
              isomorph = list(isomorph_merged_Peptides_1_2,
                              isomorph_merged_Peptides_1_3,
                              isomorph_merged_Peptides_1_4,
                              isomorph_merged_Peptides_1_5,
                              isomorph_merged_Peptides_1_6,
                              isomorph_merged_Peptides_1_7,
                              isomorph_merged_Peptides_1_8,
                              isomorph_merged_Peptides_1_9,
                              isomorph_merged_Peptides_2_3,
                              isomorph_merged_Peptides_2_4,
                              isomorph_merged_Peptides_2_5,
                              isomorph_merged_Peptides_2_6,
                              isomorph_merged_Peptides_2_7,
                              isomorph_merged_Peptides_2_8,
                              isomorph_merged_Peptides_2_9,
                              isomorph_merged_Peptides_3_4,
                              isomorph_merged_Peptides_3_5,
                              isomorph_merged_Peptides_3_6,
                              isomorph_merged_Peptides_3_7,
                              isomorph_merged_Peptides_3_8,
                              isomorph_merged_Peptides_3_9,
                              isomorph_merged_Peptides_4_5,
                              isomorph_merged_Peptides_4_6,
                              isomorph_merged_Peptides_4_7,
                              isomorph_merged_Peptides_4_8,
                              isomorph_merged_Peptides_4_9,
                              isomorph_merged_Peptides_5_6,
                              isomorph_merged_Peptides_5_7,
                              isomorph_merged_Peptides_5_8,
                              isomorph_merged_Peptides_5_9,
                              isomorph_merged_Peptides_6_7,
                              isomorph_merged_Peptides_6_8,
                              isomorph_merged_Peptides_6_9,
                              isomorph_merged_Peptides_7_8,
                              isomorph_merged_Peptides_7_9,
                              isomorph_merged_Peptides_8_9))


rownames(RES_2) <- c("proteins", "protein nodes", "peptides", "peptide nodes", "edges", "graphs",
                   "Graphs with 1 protein node",
                   "isomorphism classes",  "protein nodes", "peptide nodes", "edges",
                   "protein nodes", "peptide nodes", "edges")


RES2_2 <- RES_2
colnames(RES2_2) <- c("1_2", "1_3", "1_4", "1_5", "1_6", "1_7", "1_8", "1_9",
                      "2_3", "2_4", "2_5", "2_6", "2_7", "2_8", "2_9",
                      "3_4", "3_5", "3_6", "3_7", "3_8", "3_9",
                      "4_5", "4_6", "4_7", "4_8", "4_9",
                      "5_6", "5_7", "5_8", "5_9",
                      "6_7", "6_8", "6_9","7_8", "7_9", "8_9")
write.xlsx(RES2_2, "data/D2/D2_quant/draft_Table3_Paper.xlsx", keepNA = TRUE, row.names = TRUE)



RES3_2 <- cbind(round(rowMeans(RES_2), 0), matrixStats::rowMins(RES_2), matrixStats::rowMaxs(RES_2))
colnames(RES3_2) <- c("mean", "min", "max")





xtable(cbind(RES3, RES3_2), digits = 0)

# % latex table generated in R 4.1.0 by xtable 1.8-4 package
# % Mon Jun 28 11:52:21 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrr}
# \hline
# & mean & min & max & mean & min & max \\
# \hline
# proteins & 8110 & 7999 & 8263 & 1195 & 1168 & 1240 \\
# protein.nodes & 5393 & 5321 & 5490 & 1053 & 1024 & 1093 \\
# peptides & 18093 & 17785 & 18552 & 5578 & 5180 & 5993 \\
# peptide.nodes & 5523 & 5436 & 5638 & 1111 & 1083 & 1153 \\
# edges & 10556 & 10373 & 10800 & 1282 & 1245 & 1328 \\
# graphs & 2880 & 2833 & 2928 & 933 & 911 & 968 \\
# Graphs.with.1.protein.node & 1755 & 1720 & 1786 & 834 & 812 & 865 \\
# isomorphism.classes & 221 & 217 & 225 & 13 & 11 & 14 \\
# protein.nodes.1 & 22 & 22 & 22 & 9 & 9 & 10 \\
# peptide.nodes.1 & 18 & 17 & 20 & 9 & 8 & 10 \\
# edges.1 & 108 & 105 & 117 & 48 & 45 & 57 \\
# protein.nodes.2 & 20 & 19 & 20 & 5 & 5 & 5 \\
# peptide.nodes.2 & 25 & 24 & 27 & 10 & 10 & 11 \\
# edges.2 & 172 & 160 & 185 & 21 & 21 & 23 \\
# \hline
# \end{tabular}
# \end{table}

