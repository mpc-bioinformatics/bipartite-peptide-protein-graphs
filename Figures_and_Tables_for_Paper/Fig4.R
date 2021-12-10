library(tidyverse)
library(openxlsx)


################################################################################
#### Fig 4) barplots comparing protein nodes with and without isoforms from D3 ####

D3_fasta <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min7AA.xlsx")
D3_quant <- read.xlsx("data/D3_without_isoforms/D3_quant/table_subgraph_characteristics_D3_quant.xlsx")

D3_iso_fasta <- read.xlsx("data/D3/D3_fasta/table_subgraph_characteristics_D3_fasta_min7AA.xlsx")
D3_iso_quant <- read.xlsx("data/D3/D3_quant/table_subgraph_characteristics_D3_quant.xlsx")


### remove column containing the "comparison"
D3_quant <- D3_quant[,-1]
D3_iso_quant <- D3_iso_quant[,-1]

D3_fasta_long <- reshape2::melt(D3_fasta, id.vars = 1)
D3_quant_long <- reshape2::melt(D3_quant, id.vars = 1)
D3_iso_fasta_long <- reshape2::melt(D3_iso_fasta, id.vars = 1)
D3_iso_quant_long <- reshape2::melt(D3_iso_quant, id.vars = 1)


vars <- levels(D3_fasta_long$variable)
D_complete <- rbind(D3_fasta_long, D3_quant_long, D3_iso_fasta_long, D3_iso_quant_long)
D_complete$dataset <- rep(c("D3_fasta", "D3_quant", "D3_iso_fasta", "D3_iso_quant"),
                          times = c(nrow(D3_fasta_long), nrow(D3_quant_long), nrow(D3_iso_fasta_long), nrow(D3_iso_quant_long)))
D_complete$dataset <- factor(D_complete$dataset, levels = c("D3_fasta", "D3_iso_fasta", "D3_quant", "D3_iso_quant"))


################################################################################
### general plot settings
base_size = 5
xlab = "Data set"
ylab = "Relative frequency"
legend_titles = c("", "", "", "")
plot_margins_default = unit(c(5.5,5.5,5.5,5.5), "points")
plot_margins = plot_margins_default
legend.key.size = unit(0.3, "cm")

################################################################################
# Protein nodes

D_complete3 <- filter(D_complete, variable %in% c("Nr_prot_node_only_unique_pep",
                                                  "Nr_prot_node_uniq_and_shared_pep",
                                                  "Nr_prot_node_only_shared_pep"))
D_complete3 <- D_complete3 %>%
  dplyr::group_by(dataset, variable) %>%
  dplyr::summarise(value=sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(perc=value/sum(value))

D_complete3$variable <- droplevels(D_complete3$variable)
levels(D_complete3$variable) <- c("only unique peptides", "unique and shared peptides", "only shared peptides")


pl <- ggplot(D_complete3) +
  geom_bar(aes(x = dataset, y = perc, fill = variable), stat = "identity", position = "stack", colour = "black", size=0.3) + #
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", plot.margin = plot_margins,
        legend.key.size = legend.key.size, legend.margin=margin(t = -5),
        legend.box.margin=margin(0,0,0,0)) +
  ylab(ylab) + xlab(xlab) +
  ggtitle("Protein nodes") +
  scale_fill_manual(values = c("only unique peptides" = "grey90",
                               "unique and shared peptides" = "grey65",
                               "only shared peptides" = "grey40"),
                    name = legend_titles[[4]]) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
pl


cairo_pdf("Paper/Paper 1/figures/Figure4.pdf", height = 5/2.54, width = 4.5/2.54)
print(pl)
dev.off()




