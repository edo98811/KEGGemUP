# #' Style igraph edges based on relation_subtype
# #' @param g An igraph graph object with an edge attribute 'relation_subtype'
# #' @return The igraph graph with styled edge attributes: color, lty, arrow.mode, label
# #' @importFrom igraph edge_attr_names E
# #' @noRd
# style_edges_igraph <- function(g) {
#   stopifnot(inherits(g, "igraph"))

#   # source: https://igraph.org/r/html/1.2.7/plot.common.html
#   # edge_style_map_relation <- list(
#   #   compound = list(color = "#000000", lty = 1L, arrow.mode = 2L, label = ""),
#   #   hidden_compound = list(color = "#D3D3D3", lty = 1L, arrow.mode = 2L, label = ""),
#   #   activation = list(color = "#FF0000", lty = 1L, arrow.mode = 2L, label = ""),
#   #   inhibition = list(color = "#0000FF", lty = 1L, arrow.mode = 3L, label = ""),
#   #   expression = list(color = "#FF6B6B", lty = 1L, arrow.mode = 2L, label = ""),
#   #   repression = list(color = "#4169E1", lty = 1L, arrow.mode = 3L, label = ""),
#   #   indirect_effect = list(color = "#A9A9A9", lty = 1L, arrow.mode = 2L, label = ""),
#   #   state_change = list(color = "#808080", lty = 1L, arrow.mode = 0L, label = ""),
#   #   binding_association = list(color = "#FFB6C1", lty = 1L, arrow.mode = 0L, label = ""),
#   #   dissociation = list(color = "#696969", lty = 1L, arrow.mode = 2L, label = ""),
#   #   missing_interaction = list(color = "#FFA500", lty = 1L, arrow.mode = 2L, label = "-/-"),
#   #   phosphorylation = list(color = "#228B22", lty = 1L, arrow.mode = 2L, label = "+p"),
#   #   dephosphorylation = list(color = "#32CD32", lty = 1L, arrow.mode = 2L, label = "-p"),
#   #   glycosylation = list(color = "#9370DB", lty = 1L, arrow.mode = 2L, label = "+g"),
#   #   ubiquitination = list(color = "#FF1493", lty = 1L, arrow.mode = 2L, label = "+u"),
#   #   methylation = list(color = "#00CED1", lty = 1L, arrow.mode = 2L, label = "+m"),
#   #   others_unknown = list(color = "#DC143C", lty = 1L, arrow.mode = 2L, label = "?")
#   # )


# edge_style_map_relation <- list(
#   ECrel = list(color = "#8A2BE2", lty = 1L, arrow.mode = 2L, label = "EC"),
#   PPrel = list(color = "#FF6347", lty = 1L, arrow.mode = 2L, label = "PP"),
#   GErel = list(color = "#FFD700", lty = 1L, arrow.mode = 2L, label = "GE"),
#   PCrel = list(color = "#20B2AA", lty = 1L, arrow.mode = 2L, label = "PC"),
#   maplink = list(color = "#FF4500", lty =  5L, arrow.mode = 1L, label = "->")
# )

# edge_style_map_reaction <-
#   list(
#     reversible   = list(color = "#008000", lty = 2L, arrow.mode = 3L, label = "<->"),
#     irreversible = list(color = "#FF4500", lty = 2L, arrow.mode = 2L, label = "->")
#   )

#   group_relation_style <- list(
#     color = "transparent", lty = 2L, arrow.mode = 0L, label = ""
#   )

#   line_relation_style <- list(
#     color = "black", lty = 1L, arrow.mode = 0L, label = ""
#   )

#   n <- igraph::ecount(g)

#   # rel_sub <- tolower(igraph::E(g)$relation_subtype_name)
#   rel_sub <- igraph::E(g)$relation_type

#   for (i in seq_len(n)) {
#     if (!is.na(rel_sub[i]) && rel_sub[i] %in% names(edge_style_map_relation)) {
#       style <- edge_style_map_relation[[rel_sub[i]]]
#       igraph::E(g)$color[i] <- style$color
#       igraph::E(g)$lty[i] <- style$lty
#       igraph::E(g)$arrow.mode[i] <- style$arrow.mode
#       igraph::E(g)$label[i] <- style$label
#     }
#   }

#   reac_type <- tolower(igraph::E(g)$reaction_type)

#   for (i in seq_len(n)) {
#     if (!is.na(reac_type[i]) && reac_type[i] %in% names(edge_style_map_reaction)) {
#       style <- edge_style_map_reaction[[reac_type[i]]]
#       igraph::E(g)$color[i] <- style$color
#       igraph::E(g)$lty[i] <- style$lty
#       igraph::E(g)$arrow.mode[i] <- style$arrow.mode
#       igraph::E(g)$label[i] <- style$label
#     }
#   }

#   edge_type <- tolower(igraph::E(g)$type)

#   for (i in seq_len(n)) {
#     if (edge_type[i] == "line") {
#       igraph::E(g)$color[i] <- line_relation_style$color
#       igraph::E(g)$lty[i] <- line_relation_style$lty
#       igraph::E(g)$arrow.mode[i] <- line_relation_style$arrow.mode
#       igraph::E(g)$label[i] <- line_relation_style$label
#     } else if (edge_type[i] == "group") {
#       igraph::E(g)$color[i] <- group_relation_style$color
#       igraph::E(g)$lty[i] <- group_relation_style$lty
#       igraph::E(g)$arrow.mode[i] <- group_relation_style$arrow.mode
#       igraph::E(g)$label[i] <- group_relation_style$label;
#     }
#   }

#   g
# }


# | KEGG `type`      | Semantics                                       | igraph closest shape | visNetwork closest shape  | Rationale                                                               |
# |  | --- |  | - | --- |
# | `rectangle`      | Gene product / protein complex / ortholog group | `rectangle`          | `box`                     | Canonical node representation in both libraries                         |
# | `circle`         | Compound, glycan, small molecule                | `circle`             | `dot`                     | Circular glyph with minimal annotation                                  |
# | `roundrectangle` | Linked pathway                                  | `rectangle`          | `box` (with borderRadius) | Rounded rectangles are not native in igraph; visNetwork can approximate |
# | `line`           | Reaction, relation, or abstract connector       | `rectangle`          | `ellipse` or `dot`        | igraph has no line-node; represent as minimal node or convert to edge   |


#' Style igraph graph nodes and edges
#' @param g An igraph graph object
#' @return The styled igraph graph object
#' @noRd
style_igraph_graph <- function(g) {
  stopifnot(inherits(g, "igraph"))

  ##  Nodes
  vertices_df <- igraph::as_data_frame(g, what = "vertices")
  vertices_df <- style_vertices_igraph(vertices_df)

  # write vertex attributes back
  for (col in names(vertices_df)) {
    igraph::vertex_attr(g, col) <- vertices_df[[col]]
  }

  ##  Edges
  if (igraph::ecount(g) > 0) {
    edges_df <- igraph::as_data_frame(g, what = "edges")
    edges_df <- style_edges_igraph(edges_df)
    for (col in names(edges_df)) {
      igraph::edge_attr(g, col) <- edges_df[[col]]
    }
  }

  ##  Graph attributes
  # Preserve title and type attributes
  if (!is.null(igraph::graph_attr(g, "title"))) {
    igraph::graph_attr(g, "title") <- igraph::graph_attr(g, "title")
  }
  if (!is.null(igraph::graph_attr(g, "type"))) {
    igraph::graph_attr(g, "type") <- igraph::graph_attr(g, "type")
  }

  g
}

#' @noRd
style_edges_igraph <- function(edges_df) {
  # Define edge style maps
  edge_style_map_relation <- list(
    ECrel = list(color = "#8A2BE2", lty = 1L, arrow.mode = 2L, label = "EC"),
    PPrel = list(color = "#FF6347", lty = 1L, arrow.mode = 2L, label = "PP"),
    GErel = list(color = "#FFD700", lty = 1L, arrow.mode = 2L, label = "GE"),
    PCrel = list(color = "#20B2AA", lty = 1L, arrow.mode = 2L, label = "PC"),
    maplink = list(color = "#FF4500", lty = 5L, arrow.mode = 1L, label = "->")
  )

  edge_style_map_reaction <- list(
    reaction_substrate_reversible = list(color = "#008000", lty = 2L, arrow.mode = 1L, label = "<->"),
    reaction_product_reversible = list(color = "#008000", lty = 2L, arrow.mode = 2L, label = "<->"),
    reaction_substrate_irreversible = list(color = "#FF4500", lty = 2L, arrow.mode = 0L, label = "->"),
    reaction_product_irreversible = list(color = "#FF4500", lty = 2L, arrow.mode = 2L, label = "->")
  )

  group_relation_style <- list(
    color = "transparent", lty = 2L, arrow.mode = 0L, label = ""
  )

  line_relation_style <- list(
    color = "black", lty = 1L, arrow.mode = 0L, label = ""
  )

  # Initialize columns for styling if they don't exist
  for (col in c("color", "lty", "arrow.mode", "label")) {
    if (!(col %in% names(edges_df))) {
      edges_df[[col]] <- NA
    }
  }

  # Apply styles based on relation_type
  n <- nrow(edges_df)
  for (i in seq_len(n)) {
    rel_sub <- edges_df$relation_type[i]
    if (!is.na(rel_sub) && rel_sub %in% names(edge_style_map_relation)) {
      style <- edge_style_map_relation[[rel_sub]]
      edges_df$color[i] <- style$color
      edges_df$lty[i] <- style$lty
      edges_df$arrow.mode[i] <- style$arrow.mode
      edges_df$label[i] <- style$label
    }
  }

  # Apply styles based on reaction_type
  for (i in seq_len(n)) {
    reac_type <- edges_df$reaction_type[i]
    if (!is.na(reac_type) && reac_type %in% names(edge_style_map_reaction)) {
      style <- edge_style_map_reaction[[reac_type]]
      edges_df$color[i] <- style$color
      edges_df$lty[i] <- style$lty
      edges_df$arrow.mode[i] <- style$arrow.mode
      edges_df$label[i] <- style$label
    }
  }
  # Apply styles based on type
  for (i in seq_len(n)) {
    edge_type <- edges_df$type[i]
    if (!is.na(edge_type) && edge_type == "line") {
      edges_df$color[i] <- line_relation_style$color
      edges_df$lty[i] <- line_relation_style$lty
      edges_df$arrow.mode[i] <- line_relation_style$arrow.mode
      edges_df$label[i] <- line_relation_style$label
    } else if (!is.na(edge_type) && edge_type == "group") {
      edges_df$color[i] <- group_relation_style$color
      edges_df$lty[i] <- group_relation_style$lty
      edges_df$arrow.mode[i] <- group_relation_style$arrow.mode
      edges_df$label[i] <- group_relation_style$label
    }
  }

  return(edges_df)
}

#' Style igraph graph nodes based on KEGG graphics_type
#' @param vertices_df Data frame of graph vertices with a 'graphics_type' column
#' @return Styled vertices_df with added 'shape' and 'vertex.color' columns
#' @noRd
style_vertices_igraph <- function(vertices_df) {
  # Apply default styles based on KEGG type
  vertices_df$size <- ifelse(is.na(vertices_df$size), 25, vertices_df$size) # to check later
  vertices_df$size[vertices_df$graphics_type == "circle"] <- 10

  # Map KEGG types to shapes
  vertices_df$shape[vertices_df$graphics_type == "rectangle"] <- "vrectangle"
  vertices_df$shape[vertices_df$graphics_type == "circle"] <- "circle"
  vertices_df$shape[vertices_df$graphics_type == "roundrectangle"] <- "vrectangle"
  vertices_df$shape[vertices_df$graphics_type == "line"] <- "circle" # ellipse?
  vertices_df$shape[vertices_df$graphics_type == "ellipse"] <- "circle" # ellipse?
  vertices_df$shape[vertices_df$graphics_type == "group"] <- "circle" # ellipse?

  # Make line nodes transparent
  vertices_df$vertex.color[vertices_df$graphics_type == "line"] <- "transparent"
  vertices_df$vertex.color[vertices_df$type == "group"] <- "transparent"
  vertices_df$size[vertices_df$type == "group"] <- 2
  vertices_df$size[vertices_df$graphics_type == "line"] <- 1
  return(vertices_df)
}

#' Scale node dimensions for better visualization.
#' @param vertices_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return vertices_df with scaled x and y coordinates.
#' @noRd
scale_dimensions <- function(vertices_df, factor = 10) {
  # Scale x and y coordinates to make the graph look nicer
  vertices_df$x <- as.numeric(vertices_df$x) * factor
  vertices_df$y <- as.numeric(vertices_df$y) * factor

  return(vertices_df)
}
