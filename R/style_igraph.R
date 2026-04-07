# styling functions -------------------------------------------------------

#' Style graph nodes and edges
#'
#' Style igraph graph nodes and edges
#'
#' @details
#' The graph elements are styled keeping in mind what we require to use within
#' KEGGemUP.
#' This function is called internally by `kegg_to_graph`.
#'
#' @param g An igraph graph object
#'
#' @return The styled igraph graph object
#'
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

  return(g)
}

#' Style the edges of a graph
#'
#' Style the edges of a graph, working on the edges data frame directly
#'
#' @importFrom scales alpha
#'
#' @noRd
style_edges_igraph <- function(edges_df) {
  # Define edge style maps
  edge_style_map_relation <- list(
    ECrel = list(color = scales::alpha(.color_edge_relation_EC, 0.6),
                 lty = 2L, arrow.mode = 2L, label = "EC"),
    PPrel = list(color = scales::alpha(.color_edge_relation_PP, 0.6),
                 lty = 2L, arrow.mode = 2L, label = "PP"),
    GErel = list(color = scales::alpha(.color_edge_relation_GE, 0.6),
                 lty = 2L, arrow.mode = 2L, label = "GE"),
    PCrel = list(color = scales::alpha(.color_edge_relation_PC, 0.6),
                 lty = 2L, arrow.mode = 2L, label = "PC"),
    maplink = list(color = scales::alpha(.color_edge_relation_maplink, 0.6),
                   lty = 5L, arrow.mode = 1L, label = "maplink")
  )
  edge_style_map_reaction <- list(
    reaction_substrate_reversible = list(color = .color_edge_substrate_reversible,
                                         lty = 1L, arrow.mode = 1L, label = "\u21C4"),
    reaction_product_reversible = list(color = .color_edge_product_reversible,
                                       lty = 1L, arrow.mode = 2L, label = "\u21C4"),
    reaction_substrate_irreversible = list(color = .color_edge_substrate_irreversible,
                                           lty = 1L, arrow.mode = 0L, label = "\u2192"),
    reaction_product_irreversible = list(color = .color_edge_product_irreversible,
                                         lty = 1L, arrow.mode = 2L, label = "\u2192")
  )
  edge_style_map_type <- list(
    line  = list(color = .color_edge_line,
                 lty = 1L, arrow.mode = 0L, label = ""),
    group = list(color = .color_edge_group,
                 lty = 2L, arrow.mode = 0L, label = "")
  )

  # Apply styles based on relation_type, reaction_type and type
  edges_df <- apply_style_map(edges_df, "relation_type", edge_style_map_relation)
  edges_df <- apply_style_map(edges_df, "reaction_type", edge_style_map_reaction)
  edges_df <- apply_style_map(edges_df, "type", edge_style_map_type)

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
#'
#' @param vertices_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return vertices_df with scaled x and y coordinates.
#'
#' @noRd
scale_dimensions <- function(vertices_df, factor = 5) {

  # Validate factor input
  default_factor <- 5
  if (!is.numeric(factor) || length(factor) != 1) {
    warning("Scaling factor should be a single numeric value. Using default value of", default_factor)
    factor <- default_factor
  } else if (factor <= 0) {
    warning("Scaling factor should be positive. Using default value of", default_factor)
    factor <- default_factor
  }
  if (factor <= 0) {
    warning("Scaling factor should be positive. Using default value of", default_factor)
    factor <- default_factor
  }

  # Scale x and y coordinates to make the graph look nicer
  vertices_df$x <- as.numeric(vertices_df$x) * factor
  vertices_df$y <- as.numeric(vertices_df$y) * factor

  return(vertices_df)
}

#' Apply a style map
#'
#' @param df Data frame of edges with x and y coordinates.
#' @param column Character value
#' @param style_map Style map provided as a list
#'
#' @return Data frame of edges, but styled
#'
#' @noRd
apply_style_map <- function(df, column, style_map) {

  # What style columns are defined in the style mapping?
  style_cols <- unique(unlist(lapply(style_map, names)))
  
  # Iterate over the style map and apply styles to the data frame
  for (i in names(style_map)) {
    style <- style_map[[i]]
    indexes <- df[[column]] == i
    indexes[is.na(indexes)] <- FALSE
    for (sc in style_cols) {
      df[[sc]][indexes] <- style[[sc]]
    }
  }

  return(df)
}



#' Cleanup the title node from a KEGG igraph object
#'
#' @param g The igraph object
#'
#' @returns An igraph object, without the node originally kept as "TITLE:..."
#'
#' @noRd
#'
#' @importFrom igraph delete_vertices
#'
#' @examples
#' # TODO add a minimal one
#' ## TODO: create a faKEGG function (and graph)
cleanup_title_node <- function(g) {
  title_node <- grep(pattern = "^TITLE:", V(g)$label)
  g <- igraph::delete_vertices(g, title_node)
  return(g)
}
