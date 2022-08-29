library(tidyverse)
library(ggpubr)

source("helper_functions/duplicated_for_sparse_matrices.R")

################################################################################
#### Fig S1) barplots peptide length vs. sharedness ####



################################################################################
#### D1_fasta

load("data/D1/D1_fasta/preprocessed/01_matrix_based_on_complete_fasta_D1_fasta.RData")
inds <- duplicated.dgCMatrix(sparseM, 2)  # duplicates in columns = proteins that are collapsed into one node
sparseM_wo_duplicates <- sparseM[, !inds]
peptide_length <- nchar(rownames(sparseM_wo_duplicates))
unique <- rowSums(sparseM_wo_duplicates) == 1
D <- data.frame(peptide_length = peptide_length, unique = unique)
D <- mutate(D, status = case_when(unique ~ "unique", TRUE ~ "shared"))
D$status <- factor(D$status, levels = c("unique", "shared"))

pl1 <- ggplot(D, aes(x = peptide_length, fill = status)) +
  geom_bar() + ggtitle("D1_fasta") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + ylab("Count") + xlim(4,51) + scale_y_continuous(labels = scales::comma)
pl1


pl3 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Percentage") + ggtitle("D1_fasta") +
  theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)+
  scale_y_continuous(labels = scales::percent)
pl3_single <- pl3 + theme_bw(base_size = 25)



################################################################################
#### D1_quant

sparseM <- readRDS("data/D1/D1_quant/preprocessed/01Matrix_all_valid_peptides.rds")
inds <- base::duplicated(t(sparseM), incomparables = FALSE)
sparseM_wo_duplicates <- sparseM[, !inds]
Peptides_sequence <- gsub( " *\\(.*?\\) *", "", rownames(sparseM_wo_duplicates))
Peptides_sequence <- gsub(" ", "", Peptides_sequence, fixed = TRUE)
Peptides_sequence <- gsub(".", "", Peptides_sequence, fixed = TRUE)
ind2 <- duplicated(Peptides_sequence)
sparseM_wo_duplicates <- sparseM_wo_duplicates[!ind2, ] # collapse peptide modifications
peptide_length <- nchar(Peptides_sequence)
peptide_length[peptide_length >50] <- 50

unique <- rowSums(sparseM_wo_duplicates) == 1
D <- data.frame(peptide_length = peptide_length, unique = unique)
D <- mutate(D, status = case_when(unique ~ "unique", TRUE ~ "shared"))
D$status <- factor(D$status, levels = c("unique", "shared"))


pl2 <- ggplot(D, aes(x = peptide_length, fill = status)) +
  geom_bar() + ggtitle("D1_quant") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + ylab("Count") + xlim(4,51)
pl2



pl4 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Percentage") + ggtitle("D1_quant") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)+
  scale_y_continuous(labels = scales::percent)
pl4



################################################################################
### 2x2 Panel

pl_aligned <- cowplot::align_plots(pl1, pl2, pl3, pl4, align = "hv")
pl1_aligned <- cowplot::ggdraw(pl_aligned[[1]])
pl2_aligned <- cowplot::ggdraw(pl_aligned[[2]])
pl3_aligned <- cowplot::ggdraw(pl_aligned[[3]])
pl4_aligned <- cowplot::ggdraw(pl_aligned[[4]])

legend <- get_legend(pl1 + theme(legend.position = "bottom"))


pl_panel <- plot_grid(pl1 + theme(legend.position = "none"),
                      pl2 + theme(legend.position = "none"),
                      pl3 + theme(legend.position = "none"),
                      pl4 + theme(legend.position = "none"),
                      align = "v", labels = c("A", "B", "C", "D"), nrow = 2)
pl_panel

plot_grid(pl_panel, legend, nrow = 2, rel_heights = c(10,1))


ggsave(pl = pl_panel,"Paper/Paper 1/figures/Supplement/FigS1_600dpi.png",
       device = "png", units = "cm", height = 12.65, width = 15, dpi = 600)

