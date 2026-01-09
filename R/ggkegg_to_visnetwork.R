#' Transform a ggkegg graph to igraph or visNetwork
#'
#' @param pathway_id KEGG pathway ID (e.g., 'hsa:04110' or '04110').
#' @param return_type Output type: 'igraph' or 'visNetwork'.
#' @param scaling_factor Numeric factor to scale node sizes.
#' @param verbose Logical indicating whether to print progress messages.
#' @return An igraph or visNetwork object representing the pathway.
#' @details This function downloads the KGML file for the specified KEGG pathway,
#' then parses it to generate a graph representation using either the igraph or visNetwork package.
#' It styles nodes and edges based on their types the output can be used for
#' visualization or further analysis.
#' If differential expression results are provided,
#' they can be mapped to the nodes using the function \code{map_results_to_graph}.
#' @examples
#' pathway <- "hsa04110" # Example pathway ID
#' graph <- kegg_to_graph(pathway)
#' kegg_to_graph(pathway, return_type = "visNetwork")
#'
#' @importFrom igraph graph_from_data_frame graph_attr make_empty_graph add_vertices delete_edges E V
#'
#' @export
kegg_to_graph <- function(pathway_id, return_type = "igraph", scaling_factor = 1.5, verbose = FALSE) {
  # Check arguments
  return_type <- match.arg(return_type, choices = c("igraph", "visNetwork"), several.ok = FALSE)

  # --- 0. Validate inputs ---
  if (!is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
  }

  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)

  # --- 1. Download KGML ---
  kgml_file <- download_kgml(pathway_id, bfc_kegg)
  if (is.null(kgml_file)) {
    warning("Failed to download KGML file for pathway ID: ", pathway_id)
    return(NULL)
  }

  # --- 2. Parse nodes and edges ---
  nodes_df <- parse_kgml_entries(kgml_file)
  edges_df <- parse_kgml_edges(kgml_file)
  graph <- expand_line_nodes_and_edges(nodes_df, edges_df)
  nodes_df <- graph$nodes
  edges_df <- graph$edges

  # --- 3. Style nodes and edges ---
  nodes_df <- style_nodes(nodes_df)
  nodes_df <- add_gene_names(nodes_df)
  nodes_df <- add_compound_names(nodes_df, bfc_map)
  nodes_df <- scale_dimensions(nodes_df, factor = scaling_factor)
  nodes_df <- add_tooltip(nodes_df)
  nodes_df <- add_group(nodes_df)
  nodes_df <- nodes_df[order(nodes_df$label), ]

  if (nrow(edges_df)) {
    edges_df <- style_edges(edges_df)
    edges_df <- add_edge_tooltip(edges_df)
  }

  # --- 5. Build pathway name ---
  pathway_name <- paste0("(", pathway_id, ") ", get_pathway_name(pathway_id))

  # --- 5. Build ---
  result <- switch(return_type,
    igraph = {
      make_igraph_graph(nodes_df, edges_df, pathway_name)
    },
    visNetwork = {
      make_vis_graph(nodes_df, edges_df, pathway_name)
    },
    stop("Invalid return_type. Must be 'igraph' or 'visNetwork'.")
  )

  return(result)
}

#' Map differential expression results to nodes
#'
#' @details This functionmaps differential expression results onto the nodes of a KEGG pathway graph.
#' The pathwhay given as input must be the output of the function \code{kegg_to_graph}.
#'
#' @param g An igraph object representing the pathway.
#' @param de_results Named list of differential expression results.
#' @param return_type Output type: 'igraph' or 'visNetwork'.
#' @param feature_column Column name in de_table containing KEGG IDs
#' (if de_results is a single data.frame).
#' @param value_column Column name in de_table containing values to map
#' (if de_results is a single data.frame).
#' @param palette Color palette for node coloring (default: "RdBu").
#' @param verbose Logical indicating whether to print progress messages.
#' @return An igraph or visNetwork object with mapped results.
#' @importFrom visNetwork visIgraph visPhysics visLegend visOptions
#' @importFrom igraph as_data_frame graph_from_data_frame graph_attr permute V E
#' @details This function can be used to map the differential expression results to the graph,
#' the input of the graph must be the output of the function
#' \code{kegg_to_graph} in the igraph format. The results to be mapped can be
#' provided either as a list or as a single data.frame. If a single data.frame is provided,
#' the default column names for KEGG IDs and values are 'KEGG_ids' and 'log2FoldChange',
#' respectively, but these can be changed using the \code{feature_column} and \code{value_column} parameters.
#'
#' @examples
#' pathway <- "hsa04110" # Example pathway ID
#' graph <- kegg_to_graph(pathway, return_type = "igraph")
#' # Example differential expression results
#' de_results <- data.frame(
#'   KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
#'   log2FoldChange = c(1.5, -2.0, 0.5)
#' )
#' vis_graph <- map_results_to_graph(graph, de_results, return_type = "visNetwork")
#'
#' @export
map_results_to_graph <- function(
    g,
    de_results,
    return_type = "visNetwork",
    feature_column = NULL,
    value_column = NULL,
    palette = "RdBu",
    verbose = FALSE) {
  # Check arguments
  return_type <- match.arg(return_type, choices = c("igraph", "visNetwork"), several.ok = FALSE)

  # Check that g is an igraph object
  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object. Use 'kegg_to_graph' to create it, and don't set output type to visNetwork.")
  }

  message("Mapping differential expression results to nodes...")

  # --- 0. Validate each entry in de_results ---
  de_results <- normalize_de_results(
    de_results,
    value_column = value_column,
    feature_column = feature_column
  )

  # --- 1. Extract nodes and edges from igraph ---
  nodes_df <- igraph::as_data_frame(g, what = "vertices")
  edges_df <- igraph::as_data_frame(g, what = "edges")
  pathway_name <- igraph::graph_attr(g, "title")

  # If no valid results, return original graph with a warning
  if (is.null(de_results) || length(de_results) == 0) {
    warning("No valid differential expression results provided. Returning original graph.")
  } else {
    pathway_name <- igraph::graph_attr(g, "title")

    # --- 2. Map DE results to nodes ---
    results_combined <- combine_results_in_dataframe(de_results)
    nodes_df <- add_results_nodes(nodes_df, results_combined)

    # --- 3. Color and style nodes and edges ---
    nodes_df <- add_colors_to_nodes(nodes_df, palettes = palette)
    nodes_df <- add_tooltip(nodes_df)
  }

  # --- 4. Build ---
  result <- switch(return_type,
    igraph = {
      make_igraph_graph(nodes_df, edges_df, pathway_name)
    },
    visNetwork = {
      make_vis_graph(nodes_df, edges_df, pathway_name)
    },
    stop("Invalid return_type. Must be 'igraph' or 'visNetwork'.")
  )

  return(result)
}

#' Create a visNetwork graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return A visNetwork object representing the graph.
#' @noRd
make_vis_graph <- function(nodes_df, edges_df, pathway_name) {
  # Shapes conversion for visNetwork
  nodes_df$shape[nodes_df$shape == "vrectangle"] <- "box"
  nodes_df$shape[nodes_df$shape == "circle"] <- "dot"

  # Different handling if no edges
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = nodes_df, main = pathway_name) # if graph has no edges
  } else {
    v <- visNetwork::visNetwork(nodes = nodes_df, edges = edges_df, main = pathway_name) # if graph has edges
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)

  v <- visNetwork::visOptions(v,
    highlightNearest = list(
      enabled = FALSE,
      # degree = 2,
      hover = FALSE
    ),
    # selectedBy = "group",
    # nodesIdSelection = TRUE
  )

  v <- visNetwork::visInteraction(v,
    dragNodes = TRUE,
    multiselect = TRUE,
    selectable = TRUE
  )
  v <- visNetwork::visEvents(v,
    selectNode = "function(nodes) {
        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});
      }",
    deselectNode = "function(nodes) {
        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});
      }"
  )

  return(v)
}

#' Create an igraph graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return An igraph object representing the graph.
#' @noRd
make_igraph_graph <- function(nodes_df, edges_df, pathway_name) {
  # Shapes conversion for igraph
  nodes_df$shape[nodes_df$shape == "box"] <- "vrectangle"
  nodes_df$shape[nodes_df$shape == "dot"] <- "circle"

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

