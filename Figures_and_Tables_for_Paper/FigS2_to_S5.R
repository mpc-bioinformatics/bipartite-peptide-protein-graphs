source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

library(igraph)

################################################################################
#### Fig S2 to S5) top10 isomorphism classes ####

################################################################################
### S2a (D1_fasta)

load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS2a_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
       widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()


################################################################################
#### S2b (D1_quant)

load("data/D1/D1_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS2b_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()



################################################################################
### S3a (D2_fasta)
load("data/D2_without_isoforms/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS3a_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()




################################################################################
#### S3b (D2_quant)

load("data/D2_without_isoforms/D2_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS3b_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()



################################################################################
### S4a (D3_fasta)
load("data/D3_without_isoforms/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_without_isoforms_fasta_min7AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS4a_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()




################################################################################
#### S4b (D3_quant)

load("data/D3_without_isoforms/D3_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS4b_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()



################################################################################
### S5a (D3_fasta with isoforms)
load("data/D3/D3_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D3_fasta_min7AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS5a_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()



################################################################################
#### S5b (D3_quant with isoforms)

load("data/D3/D3_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS5b_600dpi.png"), height = 27, width = 10,
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE),
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs,
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.2,
                 legend = FALSE, vertex.size = 35, vertex.label.cex = 1.5,
                 edge.width = 3, vertex.size2 = 35, which_graphs = c(1:10),
                 margin = 0, useCanonicalPermutation = TRUE, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
     oma = c(0,0,0,0), mar = c(0,0,0,0))
legend(0.2, 0.15, legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
       xjust = 0.5, yjust = 1)
dev.off()

