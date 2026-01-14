
#' Convert KEGG pathway to graph
#' @details This function downloads the KGML file for a given KEGG pathway ID,
#' parses it, constructs an igraph object, standardizes its attributes,
#' and styles it for visualization.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param scaling_factor Numeric factor to scale node dimensions (default: 1.5).
#' @param verbose Logical indicating whether to print progress messages (default: FALSE).
#' @param simplified_graph Logical indicating whether to simplify the graph by removing KEGG-specific columns (default: TRUE).
#' @param kgml_file Optional local path to a KGML file. If provided, the function will use this file instead of downloading it.
#' @return An igraph object representing the KEGG pathway graph.
#' @export
kegg_to_graph <- function(
    pathway_id,
    scaling_factor = 1.5,
    verbose = FALSE,
    simplified_graph = TRUE,
    kgml_file = NULL) {

  # Validate pathway ID format
  if (!is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
  }

  # Setup BiocFileCache for caching downloads
  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)

  # Download KGML file if not provided
  if (is.null(kgml_file)) {
    if (verbose) {
      message("Downloading KGML file for pathway ID: ", pathway_id)
    }
    kgml_file <- download_kgml(pathway_id, bfc = bfc_kegg)
    if (is.null(kgml_file)) {
      warning("Failed to download KGML file for pathway ID: ", pathway_id)
      return(NULL)
    }
  }

  pathway_name <- paste0("(", pathway_id, ") ", get_pathway_name(pathway_id))

  g <- build_kegg_graph(kgml_file, pathway_name, bfc = bfc_map)

  g <- standardize_network(g,
    kegg_to_general_node_map(),
    kegg_to_general_edge_map(),
    node_defaults(),
    edge_defaults(),
    simplified_graph = TRUE
  )

  g <- style_igraph_graph(g, bfc_map, scaling_factor = scaling_factor)

  return(g)
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
    feature_column = NULL,
    value_column = NULL,
    palette = "RdBu",
    verbose = FALSE) {
  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }

  message("Mapping differential expression results to nodes...")

  de_results <- normalize_de_results(
    de_results,
    value_column = value_column,
    feature_column = feature_column
  )

  if (is.null(de_results) || length(de_results) == 0) {
    warning("No valid differential expression results provided. Returning original graph.")
    return(g)
  }

  # Combine results into a single data frame
  results_combined <- combine_results_in_dataframe(de_results)

  # Get current nodes
  nodes_df <- igraph::as_data_frame(g, what = "vertices")

  # Merge results into nodes
  nodes_updated <- add_results_nodes(nodes_df, results_combined)
  nodes_updated <- add_colors_to_nodes(nodes_updated, palettes = palette)

  # Order nodes to match original graph
  nodes_updated <- nodes_updated[match(igraph::V(g)$name, nodes_updated$name), ]

  # Update graph attributes
  igraph::vertex_attr(g, "de_value") <- nodes_updated$de_value[name_map[node_names]]
  igraph::vertex_attr(g, "de_source") <- nodes_updated$de_source[name_map[node_names]]
  igraph::vertex_attr(g, "color") <- nodes_updated$color[name_map[node_names]]

  return(g)
}

#' Plot KEGG pathway graph using visNetwork
#' @param g visNetwork object representing the pathway graph
#' @return visNetwork plot of the KEGG pathway
#' @importFrom igraph as_data_frame graph_attr
#' @export 
plot_visNetwork_kegg <- function(g) {
  # Convert igraph to data frames
  nodes_df <- as_data_frame(g, what = "vertices")
  edges_df <- as_data_frame(g, what = "edges")
  pathway_name <- igraph::graph_attr(g, "title")

  # Style nodes and edges
  nodes_df <- kegg_nodes_to_visNetwork(nodes_df)
  edges_df <- igraph_edges_to_visNetwork(edges_df)

  # Add tooltips
  nodes_df <- add_node_tooltip(nodes_df)
  edges_df <- add_edge_tooltip(edges_df)

  # Create visNetwork graph
  v <- make_vis_graph(nodes_df, edges_df, pathway_name)

  return(v)
}


plot_static_kegg <- function(g) {
  # ggraph(g, layout = "manual", x = x, y = y) +
  #   annotation_raster(img, xmin = 0, xmax = 50, ymin = 0, ymax = 50) +
  #   geom_edge_link() +
  #   geom_node_point()
}
