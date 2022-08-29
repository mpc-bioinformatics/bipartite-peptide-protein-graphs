library(tidyverse)
library(ggpubr)     ### for arranging plots
library(cowplot)    ### for aligning plots
library(openxlsx)


################################################################################
#### Fig 3) 2x2 Panel of barplots ####


################################################################################
# read in and preprocess data

D1_fasta <- read.xlsx("data/D1/D1_fasta/table_subgraph_characteristics_D1_fasta_min7AA.xlsx")
D1 <- read.xlsx("data/D1/D1_quant/table_subgraph_characteristics_D1_quant.xlsx")
D2_fasta <- read.xlsx("data/D2_without_isoforms/D2_fasta/table_subgraph_characteristics_D2_fasta_min6AA.xlsx")
D2 <- read.xlsx("data/D2_without_isoforms/D2_quant/table_subgraph_characteristics_D2_quant.xlsx")
D3_fasta <- read.xlsx("data/D3_without_isoforms/D3_fasta/table_subgraph_characteristics_D3_without_isoforms_fasta_min7AA.xlsx")
D3 <- read.xlsx("data/D3_without_isoforms/D3_quant/table_subgraph_characteristics_D3_quant.xlsx")


### remove column containing the "comparison"
D1 <- D1[,-1]
D2 <- D2[,-1]
D3 <- D3[,-1]

D1_fasta_long <- reshape2::melt(D1_fasta, id.vars = 1)
D1_long <- reshape2::melt(D1, id.vars = 1)
D2_fasta_long <- reshape2::melt(D2_fasta, id.vars = 1)
D2_long <- reshape2::melt(D2, id.vars = 1)
D3_fasta_long <- reshape2::melt(D3_fasta, id.vars = 1)
D3_long <- reshape2::melt(D3, id.vars = 1)

vars <- levels(D1_long$variable)
D_complete <- rbind(D1_fasta_long, D1_long, D2_fasta_long, D2_long, D3_fasta_long, D3_long)
D_complete$dataset <- rep(c("D1_fasta", "D1_quant", "D2_fasta", "D2_quant", "D3_fasta", "D3_quant"),
                          times = c(nrow(D1_fasta_long), nrow(D1_long), nrow(D2_fasta_long), nrow(D2_long), nrow(D3_fasta_long), nrow(D3_long)))
D_complete$dataset <- factor(D_complete$dataset, levels = c("D1_fasta", "D1_quant", "D2_fasta", "D2_quant", "D3_fasta", "D3_quant"))


################################################################################
### general plot settings
base_size = 5
xlab = "Data set"
ylab = "Percentage"
legend_titles = c("", "", "", "")
plot_margins_default = unit(c(5.5,5.5,5.5,5.5), "points")
plot_margins = plot_margins_default
legend.key.size = unit(0.3, "cm")

################################################################################
# Plot 1: Does graph have > 1 peptide node?

D_long_tmp <- D_complete[D_complete$variable == "has_multiple_prot",]
D_long_tmp3 <- D_long_tmp %>%
  dplyr::group_by(dataset, value) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(perc=n/sum(n))


D_long_tmp3$value <- factor(D_long_tmp3$value, labels  = c("1", "\u2265 2"))
pl1 <- ggplot(D_long_tmp3) +
  geom_bar(aes(x = dataset, y = perc, fill = value), stat = "identity", position = "stack", col = "black", size=0.3) +
  scale_fill_manual(values = c("grey90", "grey40"), name = legend_titles[[1]]) +
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", plot.margin = plot_margins,
        legend.key.size = legend.key.size, legend.margin=margin(t = -5)) +
  xlab(xlab) + ylab(ylab) +
  ggtitle("Protein nodes per graph")  +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  scale_y_continuous(labels = scales::percent)
pl1



################################################################################
# Plot 2: Does graph have >= 1 protein node\n without unique peptides?

D_long_tmp <- D_complete[D_complete$variable == "has_prot_without_unique_pep",]
D_long_tmp3 <- D_long_tmp %>%
  group_by(dataset, value) %>%
  summarise(n=n()) %>%
  group_by(dataset) %>%
  mutate(perc=n/sum(n))


D_long_tmp3$value <- factor(D_long_tmp3$value, labels  = c("0", "\u2265 1"))
pl2 <- ggplot(D_long_tmp3) +
  geom_bar(aes(x = dataset, y = perc, fill = value), stat = "identity", position = "stack", colour = "black", size=0.3) + #
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", plot.margin = plot_margins,
        legend.key.size = legend.key.size, legend.margin=margin(t = -5)) +
  ylab(ylab) + xlab(xlab) +
  ggtitle("Protein nodes without \nunique peptides per graph") +
  #ggtitle("Does graph have \u2265 1 protein \nnode without unique peptides?") +
  scale_fill_manual(values = c("grey90", "grey40"), name = legend_titles[[2]]) +  # c("none" = colours[1], "\u2265 1" = colours[3])
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))+
  scale_y_continuous(labels = scales::percent)
pl2


################################################################################
# Plot 3: Peptide nodes

D_complete2 <- filter(D_complete, variable %in% c("Nr_pep_node_unique", "Nr_pep_node_shared"))

D_complete2 <- D_complete2 %>%
  dplyr::group_by(dataset, variable) %>%
  dplyr::summarise(value=sum(value, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(perc=value/sum(value))

D_complete2$variable <- droplevels(D_complete2$variable)
levels(D_complete2$variable) <- c("unique", "shared")

pl3 <- ggplot(D_complete2) +
  geom_bar(aes(x = dataset, y = perc, fill = variable), stat = "identity", position = "stack", col = "black", size=0.3) + #
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", plot.margin = plot_margins,
        legend.key.size = legend.key.size, legend.margin=margin(t = -5)) +
  ylab(ylab) + xlab(xlab) +
  ggtitle("Uniqueness of peptide nodes") +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = legend_titles[[3]]) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), )+
  scale_y_continuous(labels = scales::percent)
pl3


################################################################################
# Plot 4: Protein nodes

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


pl4 <- ggplot(D_complete3) +
  geom_bar(aes(x = dataset, y = perc, fill = variable), stat = "identity", position = "stack", colour = "black", size=0.3) + #
  theme_bw(base_size = base_size) +
  theme(legend.position = "bottom", plot.margin = plot_margins,
        legend.key.size = legend.key.size, legend.margin=margin(t = -5),
        legend.box.margin=margin(0,0,0,0)) +
  ylab(ylab) + xlab(xlab) +
  ggtitle("Breakdown of protein nodes") +
  scale_fill_manual(values = c("only unique peptides" = "grey90",
                               "unique and shared peptides" = "grey65",
                               "only shared peptides" = "grey40"),
                    name = legend_titles[[4]]) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))+
  scale_y_continuous(labels = scales::percent)
pl4



################################################################################
### align and arrange 2x2 panel

pl_aligned <- cowplot::align_plots(pl1, pl2, pl3, pl4, align = "hv")
pl1_aligned <- cowplot::ggdraw(pl_aligned[[1]])
pl2_aligned <- cowplot::ggdraw(pl_aligned[[2]])
pl3_aligned <- cowplot::ggdraw(pl_aligned[[3]])
pl4_aligned <- cowplot::ggdraw(pl_aligned[[4]])


pl_panel <- ggpubr::ggarrange(pl1_aligned, pl2_aligned, pl3_aligned, pl4_aligned,
                              common.legend = FALSE, legend = "bottom",
                              labels = "AUTO", label.x = 0.05, label.y = 1,
                              font.label = list(size = 10, color = "black", face = "bold", family = NULL))

cairo_pdf("Paper/Paper 1/figures/Fig3.pdf", height = 10/2.54, width = 8.5/2.54)
print(pl_panel)
dev.off()

ggsave("Paper/Paper 1/figures/Fig3.tif", plot = pl_panel, width = 8.5, height = 10, device = "tiff", units = "cm", dpi = 300)
