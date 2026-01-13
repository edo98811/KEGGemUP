kegg_to_graph <- function(
    pathway_id,
    return_type = "igraph",
    scaling_factor = 1.5,
    verbose = FALSE) {

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

  pathway_name <- paste0("(", pathway_id, ") ", get_pathway_name(pathway_id))

  # --- 2. Parse nodes and edges ---
  g <- kgml_to_igraph(kgml_file, pathway_name)
  g <- standardize_igraph_graph(g, kegg_to_general_node_map, kegg_to_general_edge_map)

  # --- 3. Style nodes and edges ---
  g <- style_igraph_graph(g, bfc_map, scaling_factor = scaling_factor)

  # --- 5. Build ---
  result <- switch(return_type,
    igraph = {
      g
    },
    visNetwork = {
      g_for_vis <- convert_to_visnetwork_dfs(g)
      make_vis_graph(g_for_vis$nodes_df, g_for_vis$edges_df, pathway_name)
    }
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