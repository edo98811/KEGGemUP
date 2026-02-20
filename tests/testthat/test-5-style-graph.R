test_style_edges_igraph <- function() {

  # Create a small graph covering all styling paths
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = c(1, 2, 3, 4, 5),
      to   = c(2, 3, 4, 5, 1),

      # relation-level styles
      relation_subtype_name = c(
        "activation",        # normal relation
        "inhibition",        # normal relation
        "unknown_type",      # should remain unset (no mapping)
        "compound",          # will be overridden by reaction
        "activation"         # will be overridden by type
      ),

      # reaction-level styles
      reaction_type = c(
        NA,
        NA,
        NA,
        "reversible",        # overrides relation style
        NA
      ),

      # type-level overrides
      type = c(
        NA,
        NA,
        NA,
        NA,
        "line"               # highest priority override
      )
    ),
    vertices = data.frame(name = 1:5),
    directed = TRUE
  )

  g_styled <- style_edges_igraph(g)

  #  attribute existence 
  edge_attrs <- igraph::edge_attr_names(g_styled)
  stopifnot(all(c("color", "lty", "arrow.mode", "label") %in% edge_attrs))

  #  relation_subtype_name mapping 
  stopifnot(igraph::E(g_styled)$color[1] == "red")      # activation
  stopifnot(igraph::E(g_styled)$arrow.mode[2] == 3L)   # inhibition

  #  unknown relation subtype (no style applied) 
  stopifnot(is.na(igraph::E(g_styled)$color[3]))
  stopifnot(is.na(igraph::E(g_styled)$lty[3]))
  stopifnot(is.na(igraph::E(g_styled)$arrow.mode[3]))

  #  reaction_type overrides relation_subtype_name 
  stopifnot(igraph::E(g_styled)$lty[4] == 2L)           # reversible
  stopifnot(igraph::E(g_styled)$arrow.mode[4] == 0L)

  #  type == "line" overrides everything 
  stopifnot(igraph::E(g_styled)$color[5] == "black")
  stopifnot(igraph::E(g_styled)$lty[5] == 1L)
  stopifnot(igraph::E(g_styled)$arrow.mode[5] == 0L)
  stopifnot(igraph::E(g_styled)$label[5] == "")

  invisible(TRUE)
}


# Test style_igraph_graph
test_style_igraph_graph <- function() {
  # Create a graph with nodes having different KEGG types
  g <- igraph::graph_from_data_frame(
    data.frame(from = c(1,2), to = c(2,3), name_3 = c("expression","repression")),
    vertices = data.frame(
      name = 1:3,
      graphics_type = c("rectangle","circle","line"),
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
  v <- igraph::V(g_styled)
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
