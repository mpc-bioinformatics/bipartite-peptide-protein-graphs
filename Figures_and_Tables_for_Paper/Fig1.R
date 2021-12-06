library(igraph)
source("helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

################################################################################
#### Fig 1) Examples for bipartite graphs ####

################################################################################
### pdf version
pdf("Paper/Paper 1/figures/Figure1.pdf", width = 8.5/2.54, height = 3.44/2.54)
op <- par(mai = c(0,0,0,0))
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE),
       widths=c(3,3,2), heights=c(1,1,1))
IM1 <- matrix(c(1,0,
                1,1,
                0,1), byrow = TRUE, nrow = 3)
colnames(IM1) <- c("A", "B")
rownames(IM1) <- c("1", "2", "3")
G1 <- graph_from_incidence_matrix(IM1)
plotBipartiteGraph(G1, legend = FALSE, vertex.size = 28, vertex.label.cex = 1,
                   vertex.size2=28, edge.width = 2, three_shapes = TRUE)
text(-1,1.1, labels = "A", cex = 1.5)
text(1.1,1.1, labels = "B", cex = 1.5)
IM2 <- matrix(c(1,0,
                1,1), byrow = TRUE, nrow = 2)
colnames(IM2) <- c("A", "B")
rownames(IM2) <- c("1", "2")
G2 <- graph_from_incidence_matrix(IM2)
plotBipartiteGraph(G2, legend = FALSE, vertex.size = 28, vertex.label.cex = 1,
                   vertex.size2=28, edge.width = 2, three_shapes = TRUE)
vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.1), ylim=c(0, 0.1),
     oma = c(0,0,0,0))
legend("center", legend = c("protein", "shared peptide", "unique peptide"),
       col = vertex.color, pch = c(19, 15, 18), cex = 0.8, pt.cex = c(1.2, 1.2, 1.4))
par(op)
dev.off()
