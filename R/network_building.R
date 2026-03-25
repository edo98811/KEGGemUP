#' Build an igraph graph from a KEGG KGML file
#' @param file Path to the KGML file
#' @param pathway_name Name of the pathway
#' @param bfc_map BiocFileCache for mapping data (optional)
#' @param verbose Logical indicating whether to print verbose messages
#' @return An igraph object representing the KEGG pathway graph
#' @importFrom xml2 read_xml
#' @noRd
build_kegg_graph <- function(
  file,
  pathway_name = "Pathway",
  bfc_map = NULL,
  verbose = FALSE
) {
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
  vertices_df <- parse_kgml_nodes(
    xml,
    kegg_vertex_defaults(),
    verbose = verbose
  )

  # Parse group nodes
  vertices_df <- rbind(
    vertices_df,
    parse_kgml_groups(
      xml,
      kegg_vertex_defaults(),
      verbose = verbose
    )
  )

  # Parse line nodes
  line_nodes <-
    parse_kgml_lines(
      xml,
      kegg_vertex_defaults(),
      verbose = verbose
    )
  vertices_df <- rbind(
    vertices_df,
    line_nodes
  )

  # Parse edges
  edges_df <- parse_kgml_relations(
    xml,
    kegg_edge_defaults(),
    verbose = verbose
  )

  # Parse reactions edges
  edges_df <- rbind(
    edges_df,
    parse_kgml_reactions(
      xml,
      kegg_edge_defaults(),
      verbose = verbose
    )
  )

  # Parse line edges
  edges_df <- rbind(
    edges_df,
    parse_kgml_lines_edges(
      line_nodes,
      kegg_edge_defaults(),
      verbose = verbose
    )
  )

  # Complete reactions by adding missing substrate-product edges
  edges_df <- complete_kgml_reactions(
    vertices_df,
    edges_df,
    kegg_edge_defaults(),
    verbose = verbose
  )

  if (verbose) {
    message(
      "Total nodes: ", nrow(vertices_df), ", Total edges: ", nrow(edges_df)
    )
  }

  # I don't want to map results on line nodes
  # Indices to process
  indexes_to_map <- which(
    vertices_df$graphics_type != "line" &
      vertices_df$graphics_type != "group"
  )

  # Apply remove_kegg_prefix_str to the selected rows
  vertices_df$ids_for_mapping[indexes_to_map] <- vapply(
    vertices_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Add informations to nodes
  vertices_df <- add_node_labels(
    vertices_df,
    bfc_map,
    verbose = verbose
  )
  vertices_df <- add_reaction_labels(
    vertices_df,
    bfc_map,
    verbose = verbose
  )

  vertices_df <- add_group(vertices_df, verbose = verbose)

  if (pathway_name == "") {
    warning("Failed to retrieve pathway name; using 'Pathway' as default.")
  }

  g <- make_igraph_graph(
    vertices_df,
    edges_df,
    pathway_name,
    verbose = verbose
  )

  igraph::graph_attr(g, "type") <- "KEGG_Pathway"
  return(g)
}

#' Create an igraph graph from nodes and edges data frames
#' @param vertices_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param verbose Logical indicating whether to print verbose messages.
#' @param pathway_name Name of the pathway for the graph title.
#' @noRd
make_igraph_graph <- function(vertices_df, edges_df, pathway_name, verbose = FALSE) {
  ## TODO: rewrite
  vertices_df <- vertices_df[order(tolower(vertices_df$label), tolower(vertices_df$name)), ]
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    message("No edges in graph.")
    fake_edges <- data.frame(from = vertices_df$name[1], to = vertices_df$name[1])
    g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = vertices_df)
    g <- igraph::delete_edges(g, igraph::E(g))
  } else {
    g <- igraph::graph_from_data_frame(edges_df, directed = TRUE, vertices = vertices_df)
  }

  if (verbose) {
    message("Graph created with ", igraph::vcount(g), " vertices and ", igraph::ecount(g), " edges.")
  }

  igraph::graph_attr(g, "title") <- pathway_name

  return(g)
}

# Logic make empty graph ->
# if vertices and nodes present than use graph_from_data_frame
# If only vertices, make empty graph and add vertices and then vertices attributes


# make_igraph_graph <- function(vertices_df, edges_df, pathway_name = NULL, verbose = FALSE) {
#   # Sort vertices by label and name (case-insensitive)
#   vertices_df <- vertices_df[order(tolower(vertices_df$label), tolower(vertices_df$name)), ]

#   # Create empty graph with vertices
#   g <- igraph::make_empty_graph(n = nrow(vertices_df), directed = TRUE)
#   igraph::vertex_attr(g) <- cbind(igraph::vertex_attr(g), vertices_df)
#   igraph::set_vertex_attrs(g) <- cbind(igraph::vertex_attr(g), vertices_df)
#   igraph::V(g)$name <- vertices_df$name

#   # Add edges if available
#   if (!is.null(edges_df) && nrow(edges_df) > 0) {
#     if (nrow(edges_df) > 0) {
#       g <- igraph::add_edges(g, t(as.matrix(edges_df[, c("from", "to")])))
#       # Add edge attributes
#       for (col in setdiff(names(edges_df), c("from", "to"))) {
#         igraph::edge_attr(g, col) <- edges_df[[col]]
#       }
#         # igraph::edge_attr(g) <- cbind(igraph::edge_attr(g), edges_df)

#   } else if (verbose) {
#     warning("No edges in graph. Creating vertex-only graph.")
#   }

#   if (verbose) {
#     message("Graph created with ", igraph::vcount(g), " vertices and ", igraph::ecount(g), " edges.")
#   }

#   igraph::graph_attr(g, "title") <- pathway_name

#   return(g)
# }

# # ' Create a visNetwork graph from nodes and edges data frames
# #' @param vertices_df Data frame of nodes.
# #' @param edges_df Data frame of edges.
# #' @param pathway_name Name of the pathway for the graph title.
# #' @return A visNetwork object representing the graph.
# #' @noRd
# make_tidygraph_graph <- function(
#   vertices_df,
#   edges_df,
#   pathway_name = NULL,
#   verbose = FALSE
# ) {
#   # Sort vertices
#   vertices_df <- vertices_df[
#     order(tolower(vertices_df$label), tolower(vertices_df$name)),
#   ]

#   # Check for empty edges
#   if (nrow(edges_df) == 0 || is.null(edges_df)) {
#     warning("No edges in graph. Creating a graph with isolated nodes.")
#     edges_df <- data.frame(from = character(), to = character())
#   }

#   # Create graph
#   g <- tidygraph::tbl_graph(
#     nodes = vertices_df,
#     edges = edges_df,
#     directed = TRUE
#   )

#   # Add pathway_name as graph attribute if provided
#   if (!is.null(pathway_name)) {
#     g <- g %>%
#       activate(graph) %>%
#       mutate(name = pathway_name)
#   }

#   # Verbose output
#   if (verbose) {
#     message(
#       "Graph created with ",
#       nrow(vertices_df), " nodes and ", nrow(edges_df), " edges."
#     )
#   }

#   return(g)
# }
