# Test style_edges_igraph
test_style_edges_igraph <- function() {
  # Create a small graph with different relation subtypes
  g <- graph_from_data_frame(
    data.frame(from = c(1,2,3), to = c(2,3,1), name_3 = c("activation","inhibition","unknown_type")),
    vertices = data.frame(name = 1:3),
    directed = TRUE
  )

  g_styled <- style_edges_igraph(g)

  # Check that edge attributes were added
  edge_attrs <- edge_attr_names(g_styled)
  stopifnot(all(c("color","lty","arrow.mode","label") %in% edge_attrs))

  # Check known types
  stopifnot(E(g_styled)$color[1] == "red")       # activation
  stopifnot(E(g_styled)$arrow.mode[2] == 3)     # inhibition

  # Check unknown type mapped to others_unknown
  stopifnot(E(g_styled)$label[3] == "?")
  stopifnot(E(g_styled)$color[3] == "black")
  stopifnot(E(g_styled)$lty[3] == 2)
}

# Test style_igraph_graph
test_style_igraph_graph <- function() {
  # Create a graph with nodes having different KEGG types
  g <- graph_from_data_frame(
    data.frame(from = c(1,2), to = c(2,3), name_3 = c("expression","repression")),
    vertices = data.frame(
      name = 1:3,
      original_shape = c("rectangle","circle","line"),
      type = c("gene","compound","group"),
      size = c(NA,30,NA)
    ),
    directed = TRUE
  )

  scaling_factor <- 2
  g_styled <- style_igraph_graph(g, bfc_map = list(), scaling_factor = scaling_factor)

  # Check that node attributes were added and scaled
  node_attrs <- vertex_attr_names(g_styled)
  stopifnot(all(c("shape","color","size","label") %in% node_attrs))

  # Check node shapes
  v <- V(g_styled)
  stopifnot(v$shape[1] == "vrectangle")  # rectangle mapped to vrectangle
  stopifnot(v$shape[2] == "circle")      # circle stays circle
  stopifnot(v$shape[3] == "dot")         # line mapped to dot

  # Check node colors
  stopifnot(v$color[3] == "transparent") # line node transparent
  stopifnot(v$color[3] == "transparent") # group node transparent

  # Check node sizes
  stopifnot(v$size[1] == 25 * scaling_factor)  # NA default filled and scaled
  stopifnot(v$size[2] == 30 * scaling_factor)  # existing size scaled
  stopifnot(v$size[3] == 1 * scaling_factor)   # line node size set to 1 and scaled

  # Check edge attributes were styled
  edge_attrs <- edge_attr_names(g_styled)
  stopifnot(all(c("color","lty","arrow.mode","label") %in% edge_attrs))
}

# Test graph with no edges
test_style_graph_no_edges <- function() {
  g <- make_empty_graph(n = 2)
  g$original_shape <- c("rectangle","line")
  g$type <- c("gene","group")
  g$size <- c(NA,NA)

  g_styled <- style_igraph_graph(g, bfc_map = list(), scaling_factor = 1)

  # Node attributes added
  stopifnot(all(c("shape","color","size","label") %in% vertex_attr_names(g_styled)))

  # No edges should exist
  stopifnot(ecount(g_styled) == 0)
}
