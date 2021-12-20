#### function to plot a bipartite graph
## G: bipartite graph
## vertex.label.dist: distance of label from center of vertex
## legend: TRUE, if legend should be added
## vertex.color: vector of 3 colours for the type of vertices (proteins, shared and unique peptides)
## vertex.size/vertex.size2: size of vertices (general size or height/width of rectangles)
## vertex.label.cex: size of vertex labels
## edge.width: line width of edges
## useCanonicalPermutation: TRUE, if the canonical form of the graph is plotted (isomorph graphs will look the same)
## three_shapes: TRUE, if shared and unique peptide nodes should get different shapes

plotBipartiteGraph <- function(G, vertex.label.dist = 0, legend = TRUE,
                               vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
                               vertex.size = 15, vertex.label.cex = 1, edge.width = 1, vertex.size2=15,
                               useCanonicalPermutation = FALSE, three_shapes = FALSE, ...) {

  require(igraph)

  V(G)$type <- !V(G)$type           # switch node types so that proteins are at the top
  # 0 = proteins, 1 = peptides

  if (useCanonicalPermutation) {
    cG <- canonical_permutation(G)
    G <- permute(G, cG$labeling)
    Layout <- layout.bipartite(G)
    names_G <- character(length(V(G)))
    pos_proteins <- Layout[,1][Layout[,2] == 1]
    pos_peptides <- Layout[,1][Layout[,2] == 0]

    names_G[Layout[,2] == 1] <- LETTERS[rank(pos_proteins)]
    names_peptides <- 1:sum(Layout[,2] == 0)
    names_G[Layout[,2] == 0] <- names_peptides[rank(pos_peptides)]

    G <- set_vertex_attr(G, name = "name", value = names_G)
  }

  type <- integer(length(V(G)))
  type[!V(G)$type] <- 1                  # "protein"
  type[V(G)$type] <- 2                   # "shared peptide"
  type[V(G)$type & degree(G) == 1] <- 3  # "unique peptide"

  if (three_shapes) {
    mydiamond <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }

      symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
              stars=1.2*cbind(vertex.size, vertex.size, vertex.size, vertex.size),
              add=TRUE, inches=FALSE)
    }
    add_shape("diamond", clip= shape_noclip,
              plot=mydiamond)
    vertex.shapes = c("circle", "crectangle", "diamond")[type]
  } else {
    vertex.shapes = c("circle", "crectangle")[V(G)$type+1]
  }

  if (legend) par(mar = c(10, 4, 4, 2) + 0.1)
  plot(G, layout = layout_as_bipartite, vertex.color=vertex.color[type],
       vertex.shape = vertex.shapes,
       vertex.label.degree = c(-pi/2, pi/2)[V(G)$type+1],
       vertex.label.dist = vertex.label.dist,
       vertex.size = vertex.size, vertex.label.cex = vertex.label.cex,
       edge.width = edge.width, vertex.size2=vertex.size2, ...)

  if (legend) {
    pch = ifelse(three_shapes, c(19, 15, 18),  c(19, 15, 15))
    print(pch)
    legend(x = -1, y = -1.3, legend = c("protein", "shared peptide", "unique peptide"),
           col = vertex.color, pch = pch)

  par(mar = c(5, 4, 4, 2) + 0.1)
  }
}



#### function that converts a submatrix X into a bipartite graph
# x: element of a submatrix list
convertToBipartiteGraph <- function(x) {
  require(igraph)

  if ("list" %in% class(x)) { # class list if it contains peptide ratios
    S <- x$X
  } else   {
    S <- x
  }

  G <- graph_from_incidence_matrix(S)
  return(G)
}


#### function that calculates list of isomorph graphs
## Submatrix: Submatrix list
## matrix: TRUE, if Submatrix is list of matrices, FALSE, if it is list of graphs
calculateIsomorphList <- function(Submatrix, matrix = TRUE) {
  require(igraph)

  if (matrix) {
    Graphs <- lapply(Submatrix, convertToBipartiteGraph)
  } else {
    Graphs <- Submatrix
  }

  isomorph_list <- list()
  k <- 1

  for (i in 1:length(Graphs)) {
    print(i)
    G <- Graphs[[i]]
    if (k == 1) {isomorph_list[[k]] <- i; k <- k + 1; next}

    for (j in 1:(k-1)) {
      iso <- isomorphic(G, Graphs[[isomorph_list[[j]][1]]]) # compare graph with 1st element if each isomorphism class
      if (iso) {
        cG <- canonical_permutation(G)
        cG <- permute(G, cG$labeling)
        cG2 <- canonical_permutation(Graphs[[isomorph_list[[j]][1]]])
        cG2 <- permute(Graphs[[isomorph_list[[j]][1]]], cG2$labeling)

        iso2 <- all(V(cG)$type == V(cG2)$type)   # test if graphs are really the same (considering the node types)

        if(iso2) {
          isomorph_list[[j]] <- c(isomorph_list[[j]], i); break
        }

      }


      #  same_nr_of_prot_and_pep <- sum(V(G)$type) == sum(V(Graphs[[isomorph_list[[j]][1]]])$type) # compare number of type 1 nodes
      #  if(iso&same_nr_of_prot_and_pep){isomorph_list[[j]] <- c(isomorph_list[[j]], i); break}  # add graph to group of isomorphic graphs
    }

    if (!iso) {isomorph_list[[k]] <- i; k <- k + 1; next}                      # if it is not isomorphic to an existing class, start a new one

  }
  return(list(isomorph_list = isomorph_list, Graphs = Graphs))
}





#### function that plots list of isomorph classes, sorted by number of occurences
## isomorph_list: result of calculateIsomorphList
## Graphs: list of graphs
## path: path to save plots
## title: if TRUE, title is added with number of occurence and percentage value
## pdf: if TRUE, plot is saved in a single odf, if FALSE as multiple pngs
## cex.title: size of title
## which_graphs: ranks of classes that should be plottet (e.g. 1:10 for top 10 classes)
## mfrow: mfrow agrument of par function to arrance graphs in one plot
## save: if TRUE, graphs will be saved as pdf or png
## ...: further arguments to plotBipartiteGraph
plotIsomorphList <- function(isomorph_list, Graphs, path, title = TRUE, pdf = TRUE,
                             cex.title = 1, which_graphs = NULL, mfrow = NULL, save = TRUE,
                             title_format = "times+percent", ...) {

  ord_le_iso <- order(lengths(isomorph_list), decreasing = TRUE) # order by number of occurrences
  if (!is.null(which_graphs)) ord_le_iso <- ord_le_iso[which_graphs]
  le_iso <- lengths(isomorph_list)  # sizes of isomorph classes
  le_iso_total <- sum(le_iso)       # total number of graphs
  percentages <- round(le_iso/le_iso_total * 100, 2) # percentafes (proportion of all graphs)


  if (!is.null(mfrow)) {par(mfrow = mfrow)}
  if(pdf & save) pdf(paste0(path, ".pdf"))
  j <- 1
  for (i in ord_le_iso) {
    if (!pdf & save) png(paste0(path, "_", j, ".png"), res = 500, units = "cm", height = 15, width = 15)
    ind <- isomorph_list[[i]][1]  # plot first element for each isomorph group
    G <- Graphs[[ind]]
    types <- vertex_attr(G)$type # type = 0 peptides, type = 1 proteins
    G <- set_vertex_attr(G, name = "name", value = c(1:sum(!types), LETTERS[1:sum(types)]))

    plotBipartiteGraph(G, vertex.label.dist = 0, ...)
    if(title & title_format == "times+percent") title(paste0(le_iso[i], " times (", formatC(percentages[i], digits = 2, format = "f"), "%)"), cex.main = cex.title)
    if(title & title_format == "percent") title(paste0(formatC(percentages[i], digits = 2, format = "f"), "%"), cex.main = cex.title)
    if(!pdf & save) dev.off()
    j <- j + 1
  }
  if(pdf & save) dev.off()

  if (!is.null(mfrow)) {par(mfrow = c(1,1))}

}


