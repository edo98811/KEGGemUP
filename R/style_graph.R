#' Style igraph edges based on relation_subtype
#' @param g An igraph graph object with an edge attribute 'relation_subtype'
#' @return The igraph graph with styled edge attributes: color, lty, arrow.mode, label
#' @importFrom igraph edge_attr_names E
#' @noRd
style_edges_igraph <- function(g) {
  stopifnot("igraph" %in% class(g))
  stopifnot("name_3" %in% igraph::edge_attr_names(g))

  # Define styles
  edge_style_map <- list(
    # Default styles for relation subtypes
    compound = list(color = "black", lty = 1L, arrow.mode = 2L, label = ""),
    hidden_compound = list(color = "lightgray", lty = 1L, arrow.mode = 2L, label = ""),
    activation = list(color = "red", lty = 1L, arrow.mode = 2L, label = ""),
    inhibition = list(color = "blue", lty = 1L, arrow.mode = 3L, label = ""),
    expression = list(color = "red", lty = 2L, arrow.mode = 2L, label = ""),
    repression = list(color = "blue", lty = 2L, arrow.mode = 3L, label = ""),
    indirect_effect = list(color = "gray", lty = 2L, arrow.mode = 2L, label = ""),
    state_change = list(color = "gray", lty = 2L, arrow.mode = 0L, label = ""),
    binding_association = list(color = "black", lty = 2L, arrow.mode = 0L, label = ""),
    dissociation = list(color = "gray", lty = 2L, arrow.mode = 2L, label = ""),
    missing_interaction = list(color = "gray", lty = 2L, arrow.mode = 2L, label = "-/-"),
    phosphorylation = list(color = "black", lty = 1L, arrow.mode = 2L, label = "+p"),
    dephosphorylation = list(color = "black", lty = 1L, arrow.mode = 2L, label = "-p"),
    glycosylation = list(color = "black", lty = 1L, arrow.mode = 2L, label = "+g"),
    ubiquitination = list(color = "black", lty = 1L, arrow.mode = 2L, label = "+u"),
    methylation = list(color = "black", lty = 1L, arrow.mode = 2L, label = "+m"),
    others_unknown = list(color = "black", lty = 2L, arrow.mode = 2L, label = "?"),

    # Default style for group relations
    group_relation = list(color = "transparent", lty = 2L, arrow.mode = 0L, label = ""),

    # Default styles for reaction types
    reversible = list(color = "black", lty = 2L, arrow.mode = 0L, label = ""),
    irreversible = list(color = "black", lty = 2L, arrow.mode = 0L, label = ""),
    line = list(color = "black", lty = 1L, arrow.mode = 0L, label = "")
  )

  # Normalize relation_subtype
  # I apply this on name 3 which is reaction type or relation subtype
  rel_sub <- tolower(igraph::E(g)$name_3)
  rel_sub <- gsub("[/ ]", "_", rel_sub)
  rel_sub[is.na(rel_sub) | !(rel_sub %in% names(edge_style_map))] <- "others_unknown"

  # Vectorized assignment of edge attributes
  igraph::E(g)$color <- vapply(rel_sub, function(x) edge_style_map[[x]]$color, character(1))
  igraph::E(g)$lty <- vapply(rel_sub, function(x) edge_style_map[[x]]$lty, integer(1))
  igraph::E(g)$arrow.mode <- vapply(rel_sub, function(x) edge_style_map[[x]]$arrow.mode, integer(1))
  igraph::E(g)$label <- vapply(rel_sub, function(x) edge_style_map[[x]]$label, character(1))

  g
}


# | KEGG `type`      | Semantics                                       | igraph closest shape | visNetwork closest shape  | Rationale                                                               |
# | ---------------- | ----------------------------------------------- | -------------------- | ------------------------- | ----------------------------------------------------------------------- |
# | `rectangle`      | Gene product / protein complex / ortholog group | `rectangle`          | `box`                     | Canonical node representation in both libraries                         |
# | `circle`         | Compound, glycan, small molecule                | `circle`             | `dot`                     | Circular glyph with minimal annotation                                  |
# | `roundrectangle` | Linked pathway                                  | `rectangle`          | `box` (with borderRadius) | Rounded rectangles are not native in igraph; visNetwork can approximate |
# | `line`           | Reaction, relation, or abstract connector       | `rectangle`          | `ellipse` or `dot`        | igraph has no line-node; represent as minimal node or convert to edge   |


#' @noRd
style_igraph_graph <- function(g, bfc_map, scaling_factor = 1.5) {
  stopifnot(inherits(g, "igraph"))

  ## ---- Nodes ----
  nodes_df <- igraph::as_data_frame(g, what = "vertices")
  nodes_df <- style_nodes(nodes_df)
  nodes_df <- scale_dimensions(nodes_df, factor = scaling_factor)

  nodes_df <- nodes_df[order(nodes_df$label), , drop = FALSE]

  # write vertex attributes back
  for (col in names(nodes_df)) {
    igraph::vertex_attr(g, col) <- nodes_df[[col]]
  }

  ## ---- Edges ----
  if (igraph::ecount(g) > 0) {
    g <- style_edges_igraph(g)
  }

  g
}

style_nodes <- function(nodes_df) {
  # Apply default styles based on KEGG type
  nodes_df$size <- ifelse(is.na(nodes_df$size), 25, nodes_df$size) # to check later

  # Map KEGG types to shapes
  nodes_df$shape[nodes_df$original_shape == "rectangle"] <- "vrectangle"
  nodes_df$shape[nodes_df$original_shape == "circle"] <- "circle"
  nodes_df$shape[nodes_df$original_shape == "roundrectangle"] <- "vrectangle"
  nodes_df$shape[nodes_df$original_shape == "line"] <- "dot" # ellipse?

  # Make line nodes fully transparent
  nodes_df$color[nodes_df$original_shape == "line"] <- "transparent"
  nodes_df$color[nodes_df$type == "group"] <- "transparent"
  nodes_df$size[nodes_df$original_shape == "line"] <- 1

  return(nodes_df)
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
