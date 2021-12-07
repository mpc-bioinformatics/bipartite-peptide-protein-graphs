library(tidyverse)
library(ggpubr)

source("helper_functions/duplicated_for_sparse_matrices.R")

################################################################################
#### Fig S1 and S2) barplots peptide length vs. sharedness ####



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
  xlab("Peptide length") + ylab("Count") + xlim(4,51)
pl1
ggsave("Paper/Paper 1/figures/Supplement/FigS1a.png", device = "png",
       width = 10, height = 10, plot = pl1, dpi = 600)

pl3 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Relative Frequency") + ggtitle("D1_fasta") +
  theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)
pl3_single <- pl3 + theme_bw(base_size = 25)
ggsave("Paper/Paper 1/figures/Supplement/FigS1c.png",
       device = "png", width = 10, height = 7, plot = pl3_single, dpi = 600)


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
ggsave("Paper/Paper 1/figures/Supplement/FigS1b.png",
       device = "png", width = 10, height = 10, plot = pl2, dpi = 600)


pl4 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Relative Frequency") + ggtitle("D1_quant") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)
pl4
ggsave("Paper/Paper 1/figures/Supplement/FigS1d.png",
       device = "png", width = 10, height = 10, plot = pl4, dpi = 600)


################################################################################
### 2x2 Panel

pl_panel <- ggpubr::ggarrange(pl1, pl2, pl3, pl4,
                              common.legend = TRUE, legend = "bottom",
                              labels = "AUTO")
pl_panel


ggsave(pl = pl_panel,"Paper/Paper 1/figures/Supplement/FigS1_600dpi.png",
       device = "png", units = "cm", height = 12.65, width = 12.65, dpi = 600)



################################################################################
#### D2_fasta

load("data/D2_without_isoforms/D2_fasta/preprocessed/01_matrix_based_on_complete_fasta_D2_fasta.RData")
inds <- duplicated.dgCMatrix(sparseM, 2)
sparseM_wo_duplicates <- sparseM[, !inds]
peptide_length <- nchar(rownames(sparseM_wo_duplicates))
unique <- rowSums(sparseM_wo_duplicates) == 1
D <- data.frame(peptide_length = peptide_length, unique = unique)
D <- mutate(D, status = case_when(unique ~ "unique", TRUE ~ "shared"))
D$status <- factor(D$status, levels = c("unique", "shared"))


pl1 <- ggplot(D, aes(x = peptide_length, fill = status)) +
  geom_bar() + ggtitle("D2_fasta") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + ylab("Count") + xlim(4,51)
pl1
ggsave("Paper/Paper 1/figures/Supplement/FigS2a.png",
       device = "png", width = 10, height = 10, plot = pl1)


pl3 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Relative frequency") + ggtitle("D2_fasta") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)
pl3_single <- pl3 + theme_bw(base_size = 25)
ggsave("Paper/Paper 1/figures/Supplement/FigS2c.png",
       device = "png", width = 10, height = 7, plot = pl3_single, dpi = 600)



################################################################################
#### D2_quant

sparseM <- readRDS("data/D2_without_isoforms/D2_quant/preprocessed/01Matrix.rds")
inds <- base::duplicated(t(sparseM), incomparables = FALSE)
sparseM_wo_duplicates <- sparseM[, !inds]
peptide_length <- nchar(rownames(sparseM_wo_duplicates))
unique <- rowSums(sparseM_wo_duplicates) == 1
D <- data.frame(peptide_length = peptide_length, unique = unique)
D <- mutate(D, status = case_when(unique ~ "unique", TRUE ~ "shared"))
D$status <- factor(D$status, levels = c("unique", "shared"))


pl2 <- ggplot(D, aes(x = peptide_length, fill = status)) +
  geom_bar() + ggtitle("D2_quant") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + ylab("Count") + xlim(4,51)
pl2
ggsave("Paper/Paper 1/figures/Supplement/FigS2b.png",
       device = "png", width = 10, height = 10, plot = pl2)


pl4 <- ggplot(D, aes(x = peptide_length, fill = status, y = (..count..)/sum(..count..))) +
  geom_bar(position="fill") + ylab("Relative frequency") + ggtitle("D2_quant") + theme_bw() +
  scale_fill_manual(values = c("unique" = "grey90", "shared" = "grey40"), name = "") +
  xlab("Peptide length") + xlim(4,51)
pl4
ggsave("Paper/Paper 1/figures/Supplement/FigS2d.png",
       device = "png", width = 10, height = 10, plot = pl4)


################################################################################

pl_panel <- ggpubr::ggarrange(pl1, pl2, pl3, pl4,
                              common.legend = TRUE, legend = "bottom",
                              labels = "AUTO")
pl_panel

ggsave(pl = pl_panel,"Paper/Paper 1/figures/Supplement/FigS2_600dpi.png",
       device = "png", units = "cm", height = 12.65, width = 12.65, dpi = 600)
