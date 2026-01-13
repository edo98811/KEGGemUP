build_kegg_graph <- function(file, pathway_name = "Pathway", bfc_map = NULL) {
  xml <- read_xml(file)

  # Parse nodes (I can do that as all the functions return the same columns)
  nodes_df <- parse_kgml_nodes(xml)
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml))
  nodes_df <- rbind(nodes_df, parse_kgml_lines_nodes(xml))

  # Parse edges
  edges_df <- parse_kgml_relations(xml)
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(xml))

  # Add informations to nodes
  nodes_df <- add_compound_names(nodes_df, bfc_map)
  nodes_df <- add_gene_names(nodes_df)
  nodes_df <- add_group(nodes_df)

  if (pathway_name == "") {
    warning("Failed to retrieve pathway name; using 'Pathway' as default.")
  }
  g <- make_igraph_graph(nodes_df, edges_df, pathway_name)
}

standardize_network <- function(g, node_map, edge_map) {
  stopifnot(inherits(g, "igraph"))


  return(g)
}

#' Create an igraph graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return An igraph object representing the graph.
#' @noRd
make_igraph_graph <- function(nodes_df, edges_df, pathway_name) {
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    fake_edges <- data.frame(from = nodes_df$name[1], to = nodes_df$name[1])
    g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = nodes_df)
    g <- igraph::delete_edges(g, igraph::E(g))
  } else {
    g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
  }

  g <- igraph::permute(g, order(igraph::V(g)$label))
  igraph::graph_attr(g, "title") <- pathway_name
  return(g)
}
