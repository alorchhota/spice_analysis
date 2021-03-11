library(igraph)

coexpression_network_plot_fn = "results/coexpression.pdf"

get_5x5_net_matrix <- function() {
  n_gene = 8
  known = matrix(
    c(0, 6, 3, 4, 0, 0, 0, 0,
      6, 0, 3, 0, 0, 0, 0, 0,
      3, 3, 0, 2, 1, 0, 0, 0,
      4, 0, 2, 0, 2, 0, 0, 0,
      0, 0, 1, 2, 0, 2, 7, 0,
      0, 0, 0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 7, 0, 0, 3,
      0, 0, 0, 0, 0, 0, 3, 0),
    nrow = n_gene,
    ncol = n_gene,
    byrow = T,
    dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
                    sprintf("Gene%s", seq_len(n_gene)))
  )
  return(known)
}

net_mat = get_5x5_net_matrix()
net = igraph::graph_from_adjacency_matrix(adjmatrix = net_mat, mode = "undirected", weighted = T)
# plot(net)
V(net)$color = "deepskyblue"

pdf(coexpression_network_plot_fn)
for(i in seq_len(5)){
  plot.igraph(
    net,
    #vertex.label = V(net)$name,
    vertex.label = "",
    layout = layout.fruchterman.reingold,
    edge.color = "black",
    edge.width = E(net)$weight
  )
}
dev.off()