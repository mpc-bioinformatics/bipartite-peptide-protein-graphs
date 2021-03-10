source("helper_function/isomorph_classes_calculation_and_plotting_functions.R")

################################################################################
#### Fig S3 to S6) top10 isomorphism classes ####


################################################################################
### S3 (D1_fasta)

load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS3_600dpi.png"), height = 27, width = 12.65, 
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0)) 
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs, 
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.5, 
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
### S5 (D2_fasta)
load("data/D2/D2_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D2_fasta_min6AA_fast.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS5_600dpi.png"), height = 27, width = 12.65, 
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE), 
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0)) 
plotIsomorphList(isomorph_list = isomorph$isomorph_list, Graphs = isomorph$Graphs, 
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.5, 
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
#### S4 (D1_quant)

load("data/D1/D1_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS4_600dpi.png"), height = 27, width = 12.65, 
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE), 
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs, 
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.5, 
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
#### S6 (D2_quant)

load("data/D2/D2_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/Supplement/"

png(paste0(path, "FigS6_600dpi.png"), height = 27, width = 12.65, 
    res = 600, units = "cm")
graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE), 
                 widths=c(1,1), heights=c(1,1,1,1,1,0.7))
par(mai = c(0, 0.5, 0.5, 0), oma = c(0,0,0,0))
plotIsomorphList(isomorph_list = isomorph_all_merged_Peptides$isomorph_list, Graphs = isomorph_all_merged_Peptides$Graphs, 
                 path = path, pdf = FALSE, save = FALSE, cex.title = 1.5, 
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

