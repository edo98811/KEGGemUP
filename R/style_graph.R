#' Add gene names to gene nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @return Updated nodes data frame with gene names added to gene nodes.
#' @noRd
add_labels_kegg <- function(nodes_df) {
  # find rows that are genes (logical index)
  idx <- which(!is.na(nodes_df$type) & nodes_df$type == "gene")
  if (length(idx) == 0) {
    return(nodes_df)
  }

  # Extraction of graphic_name, handle NA
  graphics_name <- as.character(nodes_df$graphics_name)
  graphics_name[is.na(graphics_name)] <- "" # Na replaced by empty
  labels <- gsub(",.*", "", graphics_name[idx]) # take first
  labels <- trimws(labels)

  nodes_df$label[idx] <- labels
  return(nodes_df)
}

# | KEGG `type`      | Semantics                                       | igraph closest shape | visNetwork closest shape  | Rationale                                                               |
# | ---------------- | ----------------------------------------------- | -------------------- | ------------------------- | ----------------------------------------------------------------------- |
# | `rectangle`      | Gene product / protein complex / ortholog group | `rectangle`          | `box`                     | Canonical node representation in both libraries                         |
# | `circle`         | Compound, glycan, small molecule                | `circle`             | `dot`                     | Circular glyph with minimal annotation                                  |
# | `roundrectangle` | Linked pathway                                  | `rectangle`          | `box` (with borderRadius) | Rounded rectangles are not native in igraph; visNetwork can approximate |
# | `line`           | Reaction, relation, or abstract connector       | `rectangle`          | `ellipse` or `dot`        | igraph has no line-node; represent as minimal node or convert to edge   |


#' @noRd
style_igraph_graph <- function(g, bfc_map, scaling_factor = 1) {
  stopifnot(inherits(g, "igraph"))

  ## ---- Nodes ----
  nodes_df <- igraph::as_data_frame(g, what = "vertices")
  nodes_df <- style_nodes(nodes_df)
  nodes_df <- scale_dimensions(nodes_df, factor = scaling_factor)
  nodes_df <- add_tooltip(nodes_df)

  nodes_df <- nodes_df[order(nodes_df$label), , drop = FALSE]

  # write vertex attributes back
  for (col in names(nodes_df)) {
    igraph::vertex_attr(g, col) <- nodes_df[[col]]
  }

  ## ---- Edges ----
  if (igraph::ecount(g) > 0) {
    edges_df <- igraph::as_data_frame(g, what = "edges")

    edges_df <- style_edges(edges_df)
    edges_df <- add_edge_tooltip(edges_df)

    # write edge attributes back
    for (col in names(edges_df)) {
      igraph::edge_attr(g, col) <- edges_df[[col]]
    }
  }

  g
}

#' Scale node dimensions for better visualization.
#' @param nodes_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return nodes_df with scaled x and y coordinates.
#' @noRd
scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- as.numeric(nodes_df$y) * factor

  return(nodes_df)
}