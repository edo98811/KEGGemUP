#' Build an igraph graph from a KEGG KGML file
#' @param file Path to the KGML file
#' @param pathway_name Name of the pathway
#' @param bfc_map BiocFileCache for mapping data (optional)
#' @return An igraph object representing the KEGG pathway graph
#' @importFrom xml2 read_xml
#' @noRd
build_kegg_graph <- function(file, pathway_name = "Pathway", bfc_map = NULL) {
  xml <- tryCatch(
    xml2::read_xml(file),
    error = function(e) {
      stop("Failed to read XML file: ", file, "\n", conditionMessage(e))
    }
  )

  # Validate XML content
  if (is.null(xml) || length(xml2::xml_children(xml)) == 0) {
    stop("XML file is empty or invalid: ", file)
  }

  # Parse nodes
  nodes_df <- parse_kgml_nodes(xml, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml, kegg_node_defaults()))
  line_nodes <- parse_kgml_lines(xml, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, line_nodes)

  # Parse edges
  edges_df <- parse_kgml_relations(xml, kegg_edge_defaults())
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml, kegg_edge_defaults()))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults()))

  # I don't want to map results on line nodes
  # Indices to process
  indexes_to_map <- which(
    nodes_df$graphics_type != "line" & nodes_df$graphics_type != "group"
  )

  # Apply remove_kegg_prefix_str to the selected rows
  nodes_df$ids_for_mapping[indexes_to_map] <- vapply(
    nodes_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Add informations to nodes
  nodes_df <- add_node_labels(nodes_df, bfc_map)
  nodes_df <- add_reaction_labels(nodes_df, bfc_map)
  nodes_df <- add_group(nodes_df)

  if (pathway_name == "") {
    warning("Failed to retrieve pathway name; using 'Pathway' as default.")
  }

  g <- make_igraph_graph(nodes_df, edges_df)
  igraph::graph_attr(g, "title") <- pathway_name
  igraph::graph_attr(g, "type") <- "KEGG_Pathway" # open for extension with other types
  return(g)
}

#' Create an igraph graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @noRd
make_igraph_graph <- function(nodes_df, edges_df, pathway_name) {
  nodes_df <- nodes_df[order(tolower(nodes_df$label)), ]
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    fake_edges <- data.frame(from = nodes_df$name[1], to = nodes_df$name[1])
    g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = nodes_df)
    g <- igraph::delete_edges(g, igraph::E(g))
  } else {
    g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
  }

  return(g)
}
