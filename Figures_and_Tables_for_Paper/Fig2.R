source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

library(igraph)

################################################################################
#### Fig 2) top10 isomorphism classes for data set D1 ####


load("data/D1/D1_fasta/isomorph_classes/isomorph_classes_merged_Peptides_D1_fasta_min7AA_fast.RData")
load("data/D1/D1_quant/isomorph_classes/isomorph_classes_all_merged_Peptides.RData")
path <- "Paper/Paper 1/figures/"



layout_matrix <- cbind(matrix(c(1,1,2:12,12), 7, 2, byrow = TRUE), rep(13, 7), matrix(c(1, 1, 2:12,12)+13, 7, 2, byrow = TRUE))

tiff(paste0(path, "Fig2.tif"), height = 27, width = 23,  # 27, 23
     res = 300, units = "cm")
graphics::layout(layout_matrix,
                 widths=c(1,1,0.3, 1,1), heights=c(0.3,1,1,1,1,1,0.7))
par(mai = c(0.1, 0.5, 0.4, 0), oma = c(0,0,0,0),  xpd=NA)


plot(1, 1, type='n',axes=FALSE, ann=FALSE, xpd=NA)
text(x = 0.6, y = 1, labels = "A", cex = 3)

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

plot(1, 1, type='n',axes=FALSE, ann=FALSE)
abline(v = 1)


plot(1, 1, type='n',axes=FALSE, ann=FALSE, xpd=NA)
text(x = 0.6, y = 1, labels = "B", cex = 3)

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

