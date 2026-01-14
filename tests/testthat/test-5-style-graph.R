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
  edge_attrs <- igraph::edge_attr_names(g_styled)
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
  g <- igraph::graph_from_data_frame(
    data.frame(from = c(1,2), to = c(2,3), name_3 = c("expression","repression")),
    vertices = data.frame(
      name = 1:3,
      original_shape = c("rectangle","circle","line"),
      type = c("gene","compound","group"),
      size = c(NA,30,NA),
      x = c(1, 2, 3), 
      y = c(4, 5, 6),
      label = c("A", "B", "C") 
    ),
    directed = TRUE
  )

  scaling_factor <- 2
  g_styled <- style_igraph_graph(g, bfc_map = list(), scaling_factor = scaling_factor)

  # Check that node attributes were added and scaled
  node_attrs <- igraph::vertex_attr_names(g_styled)
  stopifnot(all(c("shape","color","size","label") %in% node_attrs))

  # Check node shapes
  v <- V(g_styled)
  stopifnot(v$shape[1] == "vrectangle")  # rectangle mapped to vrectangle
  stopifnot(v$shape[2] == "circle")      # circle stays circle
  stopifnot(v$shape[3] == "dot")         # line mapped to dot

  # Check node colors
  stopifnot(v$color[3] == "transparent") # line node transparent
  stopifnot(v$color[3] == "transparent") # group node transparent

  # Check node positions scaled
  stopifnot(v$x[1] == 1 * scaling_factor)
  stopifnot(v$y[1] == 4 * scaling_factor)

  # Check edge attributes were styled
  edge_attrs <- igraph::edge_attr_names(g_styled)
  stopifnot(all(c("color","lty","arrow.mode","label") %in% edge_attrs))
}

test_that("add_group correctly assigns group labels", {
  nodes <- kgml_steps$all_nodes

  nodes$label <- nodes$name

  # Run add_group
  nodes_updated <- add_group(nodes)

  # Identify group nodes
  group_nodes <- nodes_updated[nodes_updated$type == "group", ]

  # All group nodes should have non-NA group label
  expect_false(any(is.na(group_nodes$group)))

  # Components of each group node should have the same group label
  for (i in seq_len(nrow(group_nodes))) {
    group_node <- group_nodes[i, ]
    comps <- strsplit(group_node$components, ";", fixed = TRUE)[[1]]
    node_idx <- match(c(comps, group_node$name), nodes_updated$name)
    expect_true(all(nodes_updated$group[node_idx] == group_node$group))
  }

  # Check that non-NA groups have the expected value
  expect_true(all(group_nodes$group == "1;2"))
})

