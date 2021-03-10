library(ggplot2)
library(tidyverse)
library(reshape2)
library(openxlsx)

################################################################################
#### Fig 2) 4x4 panel of distribution of node types ####


################################################################################
# read in and preprocess data

D1_fasta <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min7AA.xlsx")
D1 <- read.xlsx("data/D1/D1_quant/table_subgraph_characteristics_D1_quant.xlsx")

D2_fasta <- read.xlsx("data/D2/D2_fasta/table_subgraph_characteristics_D2_fasta_min6AA.xlsx")
D2 <- read.xlsx("data/D2/D2_quant/table_subgraph_characteristics_D2_quant.xlsx")

## delete column with "comparison" for quantitative level
D1 <- D1[,-1]
D2 <- D2[,-1]

D1_fasta_long <- reshape2::melt(D1_fasta, id.vars = 1)
D1_long <- reshape2::melt(D1, id.vars = 1)
D2_fasta_long <- reshape2::melt(D2_fasta, id.vars = 1)
D2_long <- reshape2::melt(D2, id.vars = 1)

vars <- levels(D1_long$variable)
D_complete <- rbind(D1_fasta_long, D1_long, D2_fasta_long, D2_long)
D_complete$dataset <- rep(c("D1_fasta", "D1_quant", "D2_fasta", "D2_quant"), 
                          times = c(nrow(D1_fasta_long), nrow(D1_long), nrow(D2_fasta_long), nrow(D2_long)))
D_complete$dataset <- factor(D_complete$dataset, levels = c("D1_fasta", "D1_quant", "D2_fasta", "D2_quant"))

### Filter for relevant variables
D_complete2 <- D_complete[D_complete$variable %in% c("Nr_prot_node", "Nr_pep_node", "Nr_pep_node_unique", 
                                                     "Nr_pep_node_shared"), ]

outlier_limit <- 10
D_long_tmp2 <- dplyr::mutate(D_complete2, value = ifelse(value > outlier_limit, outlier_limit, value))
D_long_tmp3 <- D_long_tmp2 %>% 
  dplyr::group_by(dataset, value, variable) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::group_by(dataset, variable) %>% 
  dplyr::mutate(perc=n/sum(n))

variable.labs <- c("Protein\n nodes", "Peptide\n nodes", "Unique\n peptide\n nodes", "Shared\n peptide\n nodes")
names(variable.labs) <- levels(droplevels(D_long_tmp3$variable))

pl <- ggplot(D_long_tmp3) +
  geom_bar(aes(value, y = perc), fill = "grey", stat = "identity", col = "black", size=0.3) + 
  theme_bw(base_size = 10) + 
  ylab("Relative frequency") + xlab("Value") +
  facet_grid(variable ~ dataset, labeller = labeller(variable = variable.labs)) +
  ylim(0, 1) +
  scale_x_continuous(breaks = 0:10, labels = c(0:9, "10+")) +
  theme(axis.text.x = element_text(hjust = c(rep(0.5, 10),0.3)))
pl

ggsave("Paper/Paper 1/figures/Fig2_350dpi.jpg", plot = pl, width = 17.8, height = 8.2, device = "jpeg", dpi = 350, units = "cm")
ggsave("Paper/Paper 1/figures/Fig2_350dpi.eps", plot = pl, width = 17.8, height = 8.2, device = "eps", dpi = 350, units = "cm")
ggsave("Paper/Paper 1/figures/Fig2_1200dpi.tif", plot = pl, width = 17.8, height = 8.2, device = "tiff", dpi = 1200, units = "cm")
ggsave("Paper/Paper 1/figures/Fig2_350dpi.tif", plot = pl, width = 17.8, height = 8.2, device = "tiff", dpi = 350, units = "cm")

