#' Build an igraph graph from a KEGG KGML file
#' @param file Path to the KGML file
#' @param pathway_name Name of the pathway
#' @param bfc_map BiocFileCache for mapping data (optional)
#' @param verbose Logical indicating whether to print verbose messages
#' @return An igraph object representing the KEGG pathway graph
#' @importFrom xml2 read_xml
#' @noRd
build_kegg_graph <- function(file, pathway_name = "Pathway", bfc_map = NULL, verbose = FALSE) {
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
  vertices_df <- parse_kgml_nodes(xml, kegg_node_defaults(), verbose = verbose)
  vertices_df <- rbind(vertices_df, parse_kgml_groups(xml, kegg_node_defaults(), verbose = verbose))
  line_nodes <- parse_kgml_lines(xml, kegg_node_defaults(), verbose = verbose)
  vertices_df <- rbind(vertices_df, line_nodes)

  # Parse edges
  edges_df <- parse_kgml_relations(xml, kegg_edge_defaults(), verbose = verbose)
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml, kegg_edge_defaults(), verbose = verbose))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults(), verbose = verbose))

  if (verbose) {
    message("Total nodes: ", nrow(vertices_df), ", Total edges: ", nrow(edges_df))
  }
  # I don't want to map results on line nodes
  # Indices to process
  indexes_to_map <- which(
    vertices_df$graphics_type != "line" & vertices_df$graphics_type != "group"
  )

  # Apply remove_kegg_prefix_str to the selected rows
  vertices_df$ids_for_mapping[indexes_to_map] <- vapply(
    vertices_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Add informations to nodes
  vertices_df <- add_node_labels(vertices_df, bfc_map, verbose = verbose)
  vertices_df <- add_reaction_labels(vertices_df, bfc_map, verbose = verbose)
  vertices_df <- add_group(vertices_df, verbose = verbose)

  if (pathway_name == "") {
    warning("Failed to retrieve pathway name; using 'Pathway' as default.")
  }

  g <- make_igraph_graph(vertices_df, edges_df, verbose = verbose)
  igraph::graph_attr(g, "title") <- pathway_name
  igraph::graph_attr(g, "type") <- "KEGG_Pathway" # open for extension with other types
  return(g)
}

#' Create an igraph graph from nodes and edges data frames
#' @param vertices_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param verbose Logical indicating whether to print verbose messages.
#' @param pathway_name Name of the pathway for the graph title.
#' @noRd
make_igraph_graph <- function(vertices_df, edges_df, pathway_name, verbose = FALSE) {
  vertices_df <- vertices_df[order(tolower(vertices_df$label), tolower(vertices_df$name)), ]
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    fake_edges <- data.frame(from = vertices_df$name[1], to = vertices_df$name[1])
    g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = vertices_df)
    g <- igraph::delete_edges(g, igraph::E(g))
  } else {
    g <- igraph::graph_from_data_frame(edges_df, directed = TRUE, vertices = vertices_df)
  }

  if (verbose) {
    message("Graph created with ", igraph::vcount(g), " vertices and ", igraph::ecount(g), " edges.")
  }

  return(g)
}
