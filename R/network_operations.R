# from KGML to graph and to interactive viz -------------------------------

#' Create a graph for a KEGG pathway
#'
#' Convert KEGG pathway to igraph object
#'
#' @param pathway_id Character, KEGG pathway ID (e.g., 'hsa04110') or
#' NULL if using kgml_file
#' @param kgml_file Path to a local KGML file (optional,
#' if pathway_id is provided)
#' @param verbose Logical indicating whether to print
#' verbose messages (default: FALSE)
#'
#' @return An igraph object representing the KEGG pathway graph
#'
#' @export
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd
#'
#' @examples
#' pathway <- "hsa04110" # Example pathway ID
#' graph <- create_kegg_graph(pathway_id = pathway, verbose = TRUE)
#' plot(graph) # Plot the graph using igraph's plotting functions
create_kegg_graph <- function(pathway_id,
                              kgml_file = NULL,
                              verbose = FALSE) {
  # Validate pathway ID format
  if (is.null(kgml_file) && !is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
  }

  # Setup BiocFileCache for caching downloads
  bfc_path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(
    cache = file.path(bfc_path, "kegg_maps"),
    ask = FALSE)
  bfc_map <- BiocFileCache(
    cache = file.path(bfc_path, "mappings"),
    ask = FALSE)

  # Download KGML file if not provided
  if (is.null(kgml_file)) {
    if (verbose) {
      message("Downloading KGML file for pathway ID: ", pathway_id)
    }
    kgml_file <- retrieve_kgml(pathway_id, bfc = bfc_kegg, verbose = verbose)
    if (is.null(kgml_file)) {
      warning("Failed to download KGML file for pathway ID: ", pathway_id)
      return(NULL)
    }
  }

  # Get pathway name
  pathway_name <- paste0("(", pathway_id, ") ", get_pathway_name(pathway_id))

  # Build graph from KGML
  g <- build_kegg_graph(kgml_file, pathway_name, bfc_map = bfc_map)

  # Style graph
  g <- style_igraph_graph(g)

  # # re-sorting the vertices alphabetically
  # rank_vertices <- rank(V(g)$name)
  # g <- permute(g, rank_vertices)

  return(g)
}

#' Create an interactive visualization of KEGG pathways with visNetwork
#'
#' Plot KEGG pathway graph using visNetwork
#'
#' @details This function takes an igraph object representing a
#' KEGG pathway graph
#'  and creates a visNetwork visualization.
#' It maps node attributes to visual properties.
#'
#' @param g An igraph object representing the
#' KEGG pathway graph.
#' @param scaling_factor Numeric factor to scale
#' node sizes (default: 1.5).
#' @param relationships Character specifying which
#' relationships to include in edges
#' ("all", "reactions", "relations", "none"; default: "all").
#' @param visualization_type Character specifying the type of
#' visualization for nodes:
#'  "standard", "positions", "node_name", or "node_size" (default: "standard").
#'
#' @return A visNetwork object representing the KEGG pathway
#' graph with mapped results.
#'
#' @export
#'
#' @importFrom igraph as_data_frame graph_from_data_frame graph_attr permute V
#'
#' @examples
#' pathway <- "hsa04110" # Example pathway ID
#' graph <- create_kegg_graph(pathway_id = pathway)
#' # Example differential expression results
#' de_results <- data.frame(
#'   KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
#'   log2FoldChange = c(1.5, -2.0, 0.5)
#' )
#' graph <- map_results_to_graph(
#'   graph,
#'   de_results,
#'   feature_column = "KEGG_ids",
#'   value_column = "log2FoldChange")
#'
#' vis_graph <- render_kegg_graph(graph, scaling_factor = 1.5,
#' relationships = "all", visualization_type = "standard")
render_kegg_graph <- function(g,
                              scaling_factor = 1.5,
                              relationships = c("all", "reactions", "relations", "none"),
                              visualization_type = c("standard", "positions", "node_name", "node_size")) {

  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }

  relationships <- match.arg(relationships)
  visualization_type <- match.arg(visualization_type)
  # Convert igraph to data frames
  vertices_df <- as_data_frame(g, what = "vertices")
  edges_df <- as_data_frame(g, what = "edges")
  pathway_name <- igraph::graph_attr(g, "title")

  # Style nodes and edges
  vertices_df <- kegg_nodes_to_visNetwork(
    vertices_df,
    scaling_factor = scaling_factor,
    visualization_type = visualization_type
  )

  if (nrow(edges_df) > 0)
    edges_df <- kegg_edges_to_visNetwork(
      edges_df,
      relationships = relationships
    )

  # Add tooltips
  vertices_df <- add_node_tooltip(vertices_df)
  edges_df <- add_edge_tooltip(edges_df)

  # Create visNetwork graph
  v <- make_vis_graph(vertices_df, edges_df, pathway_name)

  return(v)
}


# mapping values onto nodes -----------------------------------------------

#' Map continuous values to graph nodes
#'
#' Map differential expression results to nodes
#'
#' @param g An igraph object representing the KEGG pathway graph.
#' @param de_results A data frame or a list of data frames containing
#'  differential expression results
#' @param feature_column Name of the column in
#' de_results that contains KEGG IDs
#' @param value_column Name of the column in de_results
#' that contains values to map
#' @param verbose Logical indicating whether to print
#' verbose messages (default: FALSE)
#' @param palette Optional color palette for mapping
#' values (default: NULL, will use a default palette)
#' @param palette_limit Optional numeric limit for the color
#'  palette (default: NULL, will be determined from data)
#' @param palettes_list Optional list of color palettes if
#' de_results is a list (default: NULL)
#' @param palettes_limits_list Optional list of numeric limits
#' for multiple palettes if de_results is a list (default: NULL)
#'
#' @return An igraph object with differential expression
#' results mapped to node attributes
#'
#' @export
#'
#' @importFrom igraph V graph_attr vertex_attr
#'
#' @details This function can be used to map the differential expression
#' results to the graph,
#' the input of the graph must be the output of the function
#' `create_kegg_graph` in the igraph format.
#' The results to be mapped can be
#' provided either as a list or as a single data.frame.
#' If a single data.frame
#' is provided, the default column names
#' that it will look for are KEGG IDs and values are
#' 'KEGG_ids' and 'log2FoldChange',
#' respectively, but these can be changed using the
#' `feature_column` and `value_column` parameters.
#'
#' @examples
#' pathway <- "hsa04110" # Example pathway ID
#' graph <- create_kegg_graph(pathway_id = pathway)
#' # Example differential expression results
#' de_results <- data.frame(
#'   KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
#'   log2FoldChange = c(1.5, -2.0, 0.5)
#' )
#' vis_graph <- map_results_to_graph(
#'   graph,
#'   de_results,
#'   feature_column = "KEGG_ids",
#'   value_column = "log2FoldChange"
#' )
map_results_to_graph <- function(g,
                                 de_results,
                                 feature_column = NULL,
                                 value_column = NULL,
                                 verbose = FALSE,
                                 palette = NULL,
                                 palette_limit = NULL,
                                 palettes_list = list(NA_character_),
                                 palettes_limits_list = c(NA_real_)) {
  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }

  if (verbose) message("Mapping differential expression results to nodes...")

  de_results <- standardize_de_results(
    de_results,
    value_column = value_column,
    feature_column = feature_column
  )

  # In case of empty or NULL results, return original graph with a warning
  if (is.null(de_results) || length(de_results) == 0) {
    warning("No valid differential expression results provided.
    Returning original graph.")
    return(g)
  }

  results_combined <- combine_results_in_dataframe(
    de_results,
    verbose = verbose
  )

  vertices_df <- igraph::as_data_frame(g, what = "vertices")

  # Merge results into nodes data frame and add colors for visualization
  nodes_updated <- add_results_nodes(
    vertices_df,
    results_combined,
    verbose = verbose
  )

  return_list <- add_colors_to_nodes(
    nodes_updated,
    palette = palette,
    palette_limit = palette_limit,
    palettes_limits_list = palettes_limits_list,
    palettes_list = palettes_list,
    verbose = verbose
  )

  nodes_updated <- return_list$vertices_df
  legend_plots <- return_list$legend_plots

  # Order nodes to match original graph
  nodes_updated <- nodes_updated[
    match(igraph::V(g)$name, nodes_updated$name), ,
    drop = FALSE
  ]

  # Update vertex attributes
  igraph::vertex_attr(g, "de_value") <- nodes_updated$de_value
  igraph::vertex_attr(g, "de_source") <- nodes_updated$de_source
  igraph::vertex_attr(g, "vertex.color") <- nodes_updated$vertex.color
  # igraph::vertex_attr(g, "text") <- nodes_updated$text
  igraph::vertex_attr(g, "de_text") <- nodes_updated$de_text
  igraph::vertex_attr(g, "de_name") <- nodes_updated$de_name
  igraph::graph_attr(g, "legend_plots") <- legend_plots

  return(g)
}


# subsetting and highlighting ---------------------------------------------

#' Subset a graph object for a KEGG pathway
#'
#' Create igraph visualization with improved layout
#'
#' @param g An igraph object to visualize.
#' Must have vertex attributes 'x' and 'y' for layout.
#' @param ids_to_include Character vector of KEGG IDs to include
#' in the subset graph.
#'
#' @details All the edges between the vertices are plotted automatically.
#'
#' @return A plot of the igraph object with improved layout.
#'
#' @export
#'
#' @importFrom igraph V induced_subgraph
#'
#' @examples
#' pathway <- "mmu00230"
#' g <- create_kegg_graph(pathway)
#' KEGG_to_include <- c("C00262", "C00385", "C00366", "C00294", "C00387",
#'                  "C01762", "C05512", "C00301", "C01185", "C00455",
#'                  "22436", "14544", "18950", "11486", "80285", "59027")
#' subg <- subset_kegg_graph(g, KEGG_to_include)
#' # plot(g)
#' # plot(subg)
subset_kegg_graph <- function(g,
                              ids_to_include) {

  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }
  if (is.null(ids_to_include) || length(ids_to_include) == 0) {
    stop("ids_to_include must be a non-empty character vector of KEGG IDs.")
  }
  # Get vertex data frame and create mapping for subgraph extraction
  vertices_df <- as_data_frame(g, what = "vertices")
  mapping <- make_mapping_df(vertices_df)

  # Subset nodes based on provided KEGG IDs
  nodes_to_include <- unique(mapping[
    mapping$matched_id %in% ids_to_include, "name"
  ])

  nodes_to_include <- nodes_to_include[!is.na(nodes_to_include)]
  if (length(nodes_to_include) == 0) {
    warning(
      "No matching nodes found for the provided KEGG IDs.",
      "Returning original graph."
    )
    return(g)
  }

  subg <- igraph::induced_subgraph(g, nodes_to_include)

  return(subg)
}

#' Highlight a subset of the KEGG graph
#'
#' Highlight a subset of the graph based on KEGG IDs
#'
#' @details This function highlights the nodes corresponding
#' to the provided KEGG IDs and fades the rest of the graph.
#' It modifies vertex attributes to achieve this effect.
#'
#' @param g An igraph object to visualize.
#' Must have vertex attributes 'x' and 'y' for layout.
#' @param ids_to_highlight Character vector of
#'  KEGG IDs to include in the highlighted subset.
#'
#' @return An igraph object with highlighted nodes
#' and faded non-highlighted nodes and edges.
#'
#' @export
#'
#' @importFrom igraph V set_vertex_attr set_edge_attr incident
#'
#' @examples
#' pathway <- "mmu00230"
#' g <- create_kegg_graph(pathway)
#' KEGG_to_include <- c("C00262", "C00385", "C00366", "C00294", "C00387",
#'                  "C01762", "C05512", "C00301", "C01185", "C00455",
#'                  "22436", "14544", "18950", "11486", "80285", "59027")
#' highlighted_subg <- highlight_kegg_graph(g, KEGG_to_include)
#'
highlight_kegg_graph <- function(g,
                                 ids_to_highlight) {
  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }

  if (is.null(ids_to_highlight) || length(ids_to_highlight) == 0) {
    stop("ids_to_highlight must be a non-empty character vector of KEGG IDs.")
  }

  # Get vertex data frame and create mapping for subgraph extraction
  vertices_df <- as_data_frame(g, what = "vertices")
  mapping <- make_mapping_df(vertices_df)

  # Subset nodes based on provided KEGG IDs
  nodes_to_highlight <- unique(
    mapping[mapping$matched_id %in% ids_to_highlight, "name"]
  )
  nodes_to_highlight <- nodes_to_highlight[!is.na(nodes_to_highlight)]

  if (length(nodes_to_highlight) == 0) {
    warning(
      "No matching nodes found for the provided KEGG IDs.",
      "Returning original graph."
    )
    return(g)
  }

  nodes_to_fade <- setdiff(igraph::V(g)$name, nodes_to_highlight)

  #  Vertex attributes
  g <- set_vertex_attr(
    g,
    "borderWidth",
    index = nodes_to_fade,
    value = 0
  )

  g <- set_vertex_attr(
    g,
    "borderWidthSelected",
    index = nodes_to_fade,
    value = 0
  )
  g <- set_vertex_attr(
    g,
    "vertex.color",
    index = nodes_to_fade,
    value = .color_node_faded
  )

  # handling a variable for the fading (maybe add the "all F" out of this and keep it in anyways?)
  g <- set_vertex_attr(
    g,
    "faded",
    value = FALSE
  )
  g <- set_vertex_attr(
    g,
    "faded",
    index = nodes_to_fade,
    value = TRUE
  )


  # handle extra the nodes with the pathway names...
  g <- set_vertex_attr(
    g,
    "font.color",
    index = nodes_to_fade,
    value = .color_pathwaynode_faded
  )

  # Edge attributes
  # Select all edges incident to any faded vertex
  edges_to_fade <- unique(
    unlist(
      lapply(
        nodes_to_fade, function(v) incident(g, v, mode = "all")
      )
    )
  )

  g <- set_edge_attr(
    g,
    "color",
    index = edges_to_fade,
    value = .color_edge_faded
  )
  g <- set_edge_attr(g, "width", index = edges_to_fade, value = 1)
  # handling a variable for the fading (maybe add the "all F" out of this and keep it in anyways?)

  ## also set the color of the edge labels to something more faded
  g <- set_edge_attr(
    g,
    "font.color",
    index = edges_to_fade,
    value = .color_edgetext_faded
  )

  return(g)
}



# exporting the graph -----------------------------------------------------

#' Export a KEGG graph
#'
#' Exporting a KEGG graph into its components, nodes and edges, as tab-separated
#' text files (having them represented as dataframes for max portability)
#'
#' @param g An igraph graph object, e.g. created with KEGGemUP
#' @param basename Character string, specifying the base name for the files to
#' write the two individual data frames, for nodes and edges
#'
#' @returns NULL, invisibly
#'
#' @importFrom igraph as_data_frame
#' @importFrom utils write.table
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#'
#' g <- create_kegg_graph(pathway_id = "hsa04110")
#' export_kegg_graph(g, basename = tempfile())
export_kegg_graph <- function(g,
                              basename) {
  # check is graph
  stopifnot(is(g, "igraph"))
  stopifnot(is.character(basename))

  # write into current directory with suffix edges and nodes
  g_nodes <- as_data_frame(g, what = "vertices")
  g_edges <- as_data_frame(g, what = "edges")

  write.table(g_nodes,
              file = paste0(basename, "_nodes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(g_edges,
              file = paste0(basename, "_edges.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("Exported graph components in ",
          basename, "_nodes.tsv and ", basename, "_edges.tsv")

  return(invisible(NULL))
}


# internal functions - from KEGG files to igraph objects ---------------------

#' Build an igraph graph from a KGML file
#'
#' Build an igraph graph from a KEGG KGML file
#'
#' @param file Character, path to the KGML file
#' @param pathway_name Name of the pathway. Character, defaults to `Pathway`
#' @param bfc_map BiocFileCache for mapping data (optional)
#' @param verbose Logical indicating whether to print verbose messages
#'
#' @return An igraph object representing the KEGG pathway graph
#'
#' @importFrom xml2 read_xml
#'
#' @noRd
build_kegg_graph <- function(file,
                             pathway_name = "Pathway",
                             bfc_map = NULL,
                             verbose = FALSE) {
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

  # strip off from the name the eventual "..." sometimes found at the end of labels
  vertices_df$graphics_name <- gsub("...", "", vertices_df$graphics_name, fixed = TRUE)

  g <- make_igraph_graph(
    vertices_df,
    edges_df,
    pathway_name,
    verbose = verbose
  )

  igraph::graph_attr(g, "type") <- "KEGG_Pathway"

  return(g)
}

#' Create an igraph object, from nodes and edges df
#'
#' Create an igraph graph from nodes and edges data frames
#'
#' @param vertices_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param verbose Logical indicating whether to print verbose messages.
#' @param pathway_name Name of the pathway for the graph title.
#'
#' @importFrom igraph E delete_edges graph_from_data_frame ecount vcount graph_attr
#'
#' @return An igraph object representing the KEGG pathway graph
#'
#' @noRd
make_igraph_graph <- function(vertices_df,
                              edges_df,
                              pathway_name,
                              verbose = FALSE) {

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

