build_kegg_graph <- function(file, pathway_name = "Pathway", bfc_map = NULL) {
  xml <- read_xml(file)

  # Parse nodes (I can do that as all the functions return the same columns)
  nodes_df <- parse_kgml_nodes(xml, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml, kegg_node_defaults()))
  line_nodes <- parse_kgml_lines(xml, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, line_nodes)

  # Parse edges
  edges_df <- parse_kgml_relations(xml, kegg_edge_defaults())
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml, kegg_edge_defaults()))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults()))

  nodes_df$KEGG_no_prefix <- vapply(nodes_df$KEGG, remove_kegg_prefix_str, character(1))

  # Add informations to nodes
  nodes_df <- add_node_labels(nodes_df, bfc_map)
  nodes_df <- add_reaction_labels(nodes_df, bfc_map)
  nodes_df <- add_group(nodes_df)

  if (pathway_name == "") {
    warning("Failed to retrieve pathway name; using 'Pathway' as default.")
  }
  g <- make_igraph_graph(nodes_df, edges_df, pathway_name)
  igraph::graph_attr(g, "title") <- pathway_name
  igraph::graph_attr(g, "type") <- pathway_name # open for extension with other types
}

#' Standardize igraph nodes and edges to default schema
#' @param g igraph object with KEGG-specific node and edge attributes
#' @param node_map named vector: KEGG attribute -> general/default attribute
#' @param edge_map named vector: KEGG attribute -> general/default attribute
#' @param node_default named list of default node attributes
#' @param edge_default named list of default edge attributes
#' @param simplified_graph logical, if TRUE keep only columns in defaults
#' @return igraph with standardized node and edge attributes
#' @noRd
standardize_network <- function(g, node_map, edge_map, node_default, edge_default, simplified_graph = TRUE) {
  stopifnot(inherits(g, "igraph"))

  ## Standardize nodes
  # Extract current node attributes
  node_attrs <- igraph::vertex_attr_names(g)
  nodes_df <- igraph::as_data_frame(g, what = "vertices")

  # Map KEGG -> general defaults
  for (kegg_attr in names(node_map)) {
    general_attr <- node_map[[kegg_attr]]
    if (kegg_attr %in% names(nodes_df)) {
      nodes_df[[general_attr]] <- nodes_df[[kegg_attr]]
    }
  }

  # Fill missing default columns
  for (col in names(node_default)) {
    if (!col %in% names(nodes_df)) {
      # This does not duplicate columns, as we only fill missing ones
      nodes_df[[col]] <- node_default[[col]]
    }
  }

  # Optionally remove KEGG-only columns
  if (simplified_graph) {
    keep_cols <- names(node_default)
    nodes_df <- nodes_df[, intersect(names(nodes_df), keep_cols), drop = FALSE]
  }

  ## Standardize edges
  edges_df <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edges_df) > 0) {
    for (kegg_attr in names(edge_map)) {
      general_attr <- edge_map[[kegg_attr]]
      if (kegg_attr %in% names(edges_df) && general_attr != "") {
        edges_df[[general_attr]] <- edges_df[[kegg_attr]]
      }
    }

    # Fill missing default columns
    for (col in names(edge_default)) {
      if (!col %in% names(edges_df)) {
        edges_df[[col]] <- edge_default[[col]]
      }
    }

    # Optionally remove KEGG-only columns
    if (simplified_graph) {
      keep_cols <- names(edge_default)
      edges_df <- edges_df[, intersect(names(edges_df), keep_cols), drop = FALSE]
    }

  } else {
    # No edges: create empty data frame with all default columns
    edges_df <- as.data.frame(lapply(edge_default, function(x) vector(mode = typeof(x), length = 0)),
      stringsAsFactors = FALSE
    )
  }

  ## Rebuild igraph with standardized attributes
  g_std <- igraph::graph_from_data_frame(
    d = edges_df,
    vertices = nodes_df,
    directed = igraph::is_directed(g)
  )

  return(g_std)
}


# #' Create an igraph graph from nodes and edges data frames
# #' @param nodes_df Data frame of nodes.
# #' @param edges_df Data frame of edges.
# #' @param pathway_name Name of the pathway for the graph title.
# #' @return An igraph object representing the graph.
# #' @noRd
# make_igraph_graph <- function(nodes_df, edges_df, pathway_name) {
#   if (nrow(edges_df) == 0 || is.null(edges_df)) {
#     warning("No edges in graph.")
#     fake_edges <- data.frame(from = nodes_df$name[1], to = nodes_df$name[1])
#     g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = nodes_df)
#     g <- igraph::delete_edges(g, igraph::E(g))
#   } else {
#     g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
#   }

#   g <- igraph::permute(g, order(igraph::V(g)$label))
#   igraph::graph_attr(g, "title") <- pathway_name
#   return(g)
# }

#' Create an igraph graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @param directed Logical, whether the graph should be directed.
#' @return An igraph object representing the graph.
#' @noRd
make_igraph_graph <- function(nodes_df, edges_df, pathway_name, directed = FALSE) {
  if (nrow(nodes_df) == 0) {
    stop("Nodes data frame is empty. Cannot create graph.")
  }

  # Sort nodes by label for consistent vertex order
  nodes_df <- nodes_df[order(nodes_df$label), , drop = FALSE]

  if (is.null(edges_df) || nrow(edges_df) == 0) {
    # Create empty graph and add vertices
    g <- igraph::graph.empty(n = nrow(nodes_df), directed = directed)
    igraph::vertex_attr(g) <- as.list(nodes_df)
  } else {
    g <- igraph::graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = directed)
  }

  # Set pathway title
  igraph::graph_attr(g, "title") <- pathway_name
  return(g)
}
