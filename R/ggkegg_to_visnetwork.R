#' Transform a ggkegg graph to igraph or visNetwork
#'
#' @param pathway_id KEGG pathway ID (e.g., 'hsa:04110' or '04110').
#' @param return_type Output type: 'igraph' or 'visNetwork'.
#' @param scaling_factor Numeric factor to scale node sizes.
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
kegg_to_graph <- function(pathway_id, return_type = "igraph", scaling_factor = 1.5) {
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
#' @return An igraph or visNetwork object with mapped results.
#' @importFrom visNetwork visIgraph visPhysics visLegend visOptions
#' @importFrom igraph as_data_frame graph_from_data_frame graph_attr permute V E
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
    palette = "RdBu") {
  # Check arguments
  return_type <- match.arg(return_type, choices = c("igraph", "visNetwork"), several.ok = FALSE)

  # Check that g is an igraph object
  if (!inherits(g, "igraph")) {
    stop("Input graph 'g' must be an igraph object.")
  }

  message("Mapping differential expression results to nodes...")

  # --- 0. Validate each entry in de_results ---
  if (!is.null(de_results)) {
    # If input is a data.frame, convert to default named list
    if (inherits(de_results, "data.frame")) {
      de_results <- list(de_input = list(de_table = de_results, value_column = ifelse(is.null(value_column),
        "log2FoldChange", value_column
      ), feature_column = ifelse(is.null(feature_column),
        "KEGG_ids", feature_column
      )))
    }

    # Check that de_results is a named list
    if (!is.list(de_results) || is.null(names(de_results)) || any(names(de_results) ==
      "")) {
      warning("de_results must be a named list or NULL. Ignoring de_results.")
      de_results <- NULL
    }

    # Keep only valid entries
    de_results <- de_results[vapply(names(de_results), function(name) {
      is_valid_de_entry(
        de_results[[name]],
        name
      )
    }, logical(1))]
  }

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
  nodes_df$shape[nodes_df$shape == "vrectangle"] <-  "box"
  nodes_df$shape[nodes_df$shape == "circle"] <-  "dot"

  # Different handling if no edges
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = nodes_df, main = pathway_name) # if graph has no edges
  } else {
    v <- visNetwork::visNetwork(nodes = nodes_df, edges = edges_df, main = pathway_name) # if graph has edges
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)

  v <- visNetwork::visOptions(v, highlightNearest = list(
    enabled = TRUE, degree = 2,
    hover = TRUE
  ), selectedBy = "group", nodesIdSelection = TRUE)

  v <- visNetwork::visInteraction(v, dragNodes = TRUE)

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
  nodes_df$shape[nodes_df$shape == "box"] <-  "vrectangle"
  nodes_df$shape[nodes_df$shape == "dot"] <-  "circle"

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

#' Add group information to nodes based on 'undefined' groups.
#' @param nodes_df Data frame of nodes with columns: id, kegg_name, components.
#' @return nodes_df with updated 'group' column.
#' @noRd
add_group <- function(nodes_df) {
  # The nodes that have as kegg name 'undefined' are group nodes
  undefined_idx <- which(!is.na(nodes_df$kegg_name) & nodes_df$kegg_name == "undefined")

  # If there are no undefined nodes, return original df
  if (length(undefined_idx) == 0) {
    return(nodes_df)
  }

  # Subset undefined nodes
  undefined_nodes <- nodes_df[undefined_idx, , drop = FALSE]

  # Iterate over undefined nodes using indeces
  for (i in seq_len(nrow(undefined_nodes))) {
    # If no components in group (empty), skip
    if (is.na(undefined_nodes$components[i]) || undefined_nodes$components[i] ==
      "") {
      next
    } # If group is NA

    # Get component ids and add the group label to them, add the group node
    # itself to this list
    ids <- strsplit(undefined_nodes$components[i], ";", fixed = TRUE)[[1]]
    ids <- append(ids, undefined_nodes$id[i])

    # Make group label group_label <- paste0('group_',
    # undefined_nodes$id[i])
    
    group_elements <- nodes_df$label[nodes_df$id %in% ids]
    group_label <- paste(group_elements[1:length(group_elements)-1], collapse = ";")

    # Assign group label to nodes_df
    nodes_df[nodes_df$id %in% ids, "group"] <- group_label
  }

  return(nodes_df)
}

#' Scale node dimensions for better visualization.
#' @param nodes_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return nodes_df with scaled x and y coordinates.
#' @noRd
scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- -as.numeric(nodes_df$y) * factor # Invert y-axis

  return(nodes_df)
}

#' Add tooltips to nodes for visNetwork visualization.
#' @param nodes_df Data frame of nodes with columns: KEGG, label, source, value.
#' @return nodes_df with added 'title' column for tooltips.
#' @details The tooltip includes a button to the specific KEGG entry page.
#' If multiple KEGG IDs are present, they are concatenated with '+' in the URL.
#' It also adds information about the node name, source of differential
#' expression data, and value.
#' @noRd
add_tooltip <- function(nodes_df) {

  button_html <- ifelse(
    is.na(nodes_df$kegg_name) | is.na(nodes_df$link) | nodes_df$kegg_name == "",
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", nodes_df$link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#337ab7; border-color:#2e6da4;'>",
      "KEGG entry",
      "</button></a></div>"
    )
  )
  nodes_df$title <- ifelse(
    nodes_df$kegg_name == "undefined",

    # Group placeholder node
    paste0(
      "<table>",
      "<tr><th align='left'>Group</th><td>",
      ifelse(
        is.na(nodes_df$group) | nodes_df$group == "",
        "Not part of any group",
        nodes_df$group
      ),
      "</td></tr>",
      "</table>"
    ),

    # Regular node
    paste0(
      "<table>",
      "<tr><th align='left'>KEGG Name</th><td>",
      ifelse(
        nchar(nodes_df$kegg_name) > 50,
        substr(nodes_df$kegg_name, 1, 50),
        nodes_df$kegg_name
      ),
      "</td></tr>",
      "<tr><th align='left'>Name</th><td>",
      ifelse(is.na(nodes_df$graphics_name), "", nodes_df$graphics_name),
      "</td></tr>",
      "<tr><th align='left'>Source</th><td>",
      ifelse(is.na(nodes_df$source), "", nodes_df$source),
      "</td></tr>",
      "<tr><th align='left'>Value</th><td>",
      ifelse(
        is.na(nodes_df$plot_value), "",
        format(round(as.numeric(nodes_df$plot_value), 3), nsmall = 3)
      ),
      "</td></tr>",
      "<tr><th align='left'>Group</th><td>",
      ifelse(
        is.na(nodes_df$group) | nodes_df$group == "",
        "Not belonging to any group",
        nodes_df$group
      ),
      "</td></tr>",
      "</table>",
      button_html
    )
  )


  return(nodes_df)
}

#' Style nodes based on their type for visNetwork visualization.
#' @param nodes_df Data frame of nodes with a column 'type'.
#' @param node_size_multiplier Numeric factor to scale node sizes (default: 1.2).
#' @return nodes_df with added visual styling columns:
#' shape, fixed, widthConstraint, heightConstraint, size.
#' @noRd
style_nodes <- function(nodes_df, node_size_multiplier = 1.2) {
  # Base visual settings (for vinetwork)
  nodes_df$shape[nodes_df$type == "compound"] <-  "dot"
  nodes_df$shape[nodes_df$type != "compound"] <-  "box"
  # Set size constraints for non-compound nodes (compute numeric vectors
  # first)
  widths_num <- as.numeric(nodes_df$width) * node_size_multiplier
  heights_num <- as.numeric(nodes_df$height) * node_size_multiplier
  non_comp_idx <- which(!is.na(nodes_df$type) & nodes_df$type != "compound")

  if (length(non_comp_idx) > 0) {
    nodes_df$widthConstraint[non_comp_idx] <- widths_num[non_comp_idx]
    nodes_df$heightConstraint[non_comp_idx] <- heights_num[non_comp_idx]
  }

  # Group nodes set dimension to one (very small)
  undef_idx <- which(!is.na(nodes_df$kegg_name) & nodes_df$kegg_name == "undefined")
  if (length(undef_idx) > 0) {
    nodes_df$widthConstraint[undef_idx] <- 1
    nodes_df$heightConstraint[undef_idx] <- 1
  }
  dot_idx <- which(nodes_df$shape == "dot")
  if (length(dot_idx) > 0) {
    nodes_df$size[dot_idx] <- 7
  }
  return(nodes_df)
}


#' Style edges based on their relation_subtype for visNetwork visualization.
#' @param edges_df Data frame of edges with a column 'relation_subtype'.
#' @return edges_df with added visual styling columns: color, dashes, arrows, label.
#' @noRd
style_edges <- function(edges_df) {
  # possible relation_subtypes and their styles
  # name	value	ECrel	PPrel	GErel	Explanation
  # compound	Entry element id attribute value for compound.	*	*		shared with two successive reactions (ECrel) or intermediate of two interacting proteins (PPrel)
  # hidden compound	Entry element id attribute value for hidden compound.	*			shared with two successive reactions but not displayed in the pathway map
  # activation	-->		*		positive and negative effects which may be associated with molecular information below
  # inhibition	--|		*
  # expression	-->			*	interactions via DNA binding
  # repression	--|			*
  # indirect effect	..>		*	*	indirect effect without molecular details
  # state change	...		*		state transition
  # binding/association	---		*		association and dissociation
  # dissociation	-+-		*
  # missing interaction	-/-		*	*	missing interaction due to mutation, etc.
  # phosphorylation	+p		*		molecular events
  # dephosphorylation	-p		*
  # glycosylation	+g		*
  # ubiquitination	+u		*
  # methylation	+m		*

  edge_style_map <- list(
    compound = list(color = "black", dashes = FALSE, arrows = "to", label = ""),
    hidden_compound = list(color = "lightgray", dashes = FALSE, arrows = "to", label = ""),
    activation = list(color = "red", dashes = FALSE, arrows = "to", label = ""),
    inhibition = list(color = "blue", dashes = FALSE, arrows = "tee", label = ""),
    expression = list(color = "red", dashes = TRUE, arrows = "to", label = ""),
    repression = list(color = "blue", dashes = TRUE, arrows = "tee", label = ""),
    indirect_effect = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    state_change = list(color = "gray", dashes = TRUE, arrows = "", label = ""),
    binding_association = list(color = "black", dashes = TRUE, arrows = "", label = ""),
    dissociation = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    missing_interaction = list(color = "gray", dashes = TRUE, arrows = "to", label = "-/-"),
    phosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+p"),
    dephosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "-p"),
    glycosylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+g"),
    ubiquitination = list(color = "black", dashes = FALSE, arrows = "to", label = "+u"),
    methylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+m"),
    others_unknown = list(color = "black", dashes = TRUE, arrows = "to", label = "?"),
    group_relation = list(color = "transparent", dashes = TRUE, arrows = "", label = "")
  )

  # https://builtin.com/data-science/and-in-r#:~:text=The%20single%20sign%20version%20%7C%20returns,first%20element%20of%20each%20vector.
  edges_df$relation_subtype <- tolower(edges_df$relation_subtype)
  edges_df$relation_subtype <- gsub("[/ ]", "_", edges_df$relation_subtype)
  edges_df$relation_subtype[is.na(edges_df$relation_subtype) |
    !(edges_df$relation_subtype %in% names(edge_style_map))] <- "others_unknown"

  # Vectorized assignment
  edges_df$color <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$color, character(1))
  edges_df$dashes <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$dashes, logical(1))
  edges_df$arrows <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$arrows, character(1))
  edges_df$label <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$label, character(1))


  return(edges_df)
}

#' Add tooltips to edges for visNetwork visualization.
#' @param edges_df Data frame of edges with columns: relation_subtype, type, label.
#' @return edges_df with added 'title' column for tooltips.
#' @noRd
add_edge_tooltip <- function(edges_df) {
  edges_df$title <- paste0(
    "relation_subtype: ", edges_df$relation_subtype, "<br>",
    "Type: ", edges_df$relation_type, "<br>",
    "Label: ", ifelse(edges_df$label == "", "N/A", edges_df$label)
  )
  return(edges_df)
}


#' Add results from combined results data frame to nodes data frame
#' @param nodes_df Data frame of nodes with a column 'KEGG'
#' @param results_combined Data frame with combined results containing columns:
#' KEGG, value, source
#' @return Updated nodes data frame with added columns: value, color, source, text
#' @noRd
add_results_nodes <- function(nodes_df, results_combined) {
  # If results is empty then return the original df
  if (is.null(results_combined)) {
    return(nodes_df)
  }

  warn <- FALSE
  # Node mapping has this structure:
  # id   | KEGG (id is the node id)
  # 1      1111
  # 2      2222
  # 2      3333
  # 3    K00001
  # 4  tst00002
  # 5 undefined
  # 6      1111

  # results_combined has this structure:
  # KEGG     | plot_value | source
  # 1111        2.5        de1
  # 2222       -1.2        de2
  # 3333        0.5        de1

  # Ensure required columns exist
  required_nodes <- c("id", "KEGG", "plot_value", "source", "text")
  required_results <- c("KEGG", "plot_value", "source")

  stopifnot(
    all(required_nodes %in% names(nodes_df)),
    all(required_results %in% names(results_combined))
  )

  # Explode KEGG IDs per node -------------------------------------------
  mapping <- do.call(
    rbind,
    lapply(seq_len(nrow(nodes_df)), function(i) {
      data.frame(
        id = nodes_df$id[i],
        KEGG = strsplit(nodes_df$KEGG[i], ";", fixed = TRUE)[[1]],
        stringsAsFactors = FALSE
      )
    })
  )

  # Join results onto exploded mapping ----------------------------------
  mapping <- merge(
    mapping,
    results_combined,
    by = "KEGG",
    all.x = FALSE,
    sort = FALSE
  )

  if (nrow(mapping) == 0) {
    return(nodes_df)
  }

  # Detect multiple matches per node ------------------------------------
  match_counts <- table(mapping$id)
  warn_nodes <- names(match_counts[match_counts > 1])
  warn <- length(warn_nodes) > 0

  # Assign first match only -----------------------------
  first_hits <- mapping[!duplicated(mapping$id), ]
  idx <- match(first_hits$id, nodes_df$id)

  nodes_df$plot_value[idx] <- ifelse(
    is.na(nodes_df$plot_value[idx]),
    first_hits$plot_value,
    nodes_df$plot_value[idx]
  )

  nodes_df$source[idx] <- ifelse(
    is.na(nodes_df$source[idx]),
    first_hits$source,
    nodes_df$source[idx]
  )

  # Append text for all matches -----------------------------------------
  sep <- ","
  mapping$text_append <- paste0(
    "Source: ", mapping$source,
    sep, "Value: ", mapping$plot_value,
    sep, "Id: ", mapping$KEGG,
    ";"
  )

  text_by_node <- tapply(
    mapping$text_append,
    mapping$id,
    paste0,
    collapse = ""
  )

  text_idx <- match(names(text_by_node), nodes_df$id)
  nodes_df$text[text_idx] <- paste0(
    nodes_df$text[text_idx],
    text_by_node
  )

  # --- 6. Warn if necessary ----------------------------------------------------
  if (warn) {
    warning(
      "Some nodes had multiple matching KEGG IDs; ",
      "only the first match was used for plot_value/source."
    )
  }

  return(nodes_df)
}

#' Combine multiple differential expression results into a single data frame
#' @param results_list A named list where each element is a differential
#' expression result containing a `data.frame`
#' (de_table), value column name (value_column), and feature column name (feature_column)
#' @return A combined data frame with columns: KEGG, value, source
#' @noRd
combine_results_in_dataframe <- function(results_list) {
  # If no results provided, return NULL
  if (is.null(results_list) || length(results_list) == 0) {
    return(NULL)
  }

  # Combine all results into a single data frame
  results <- lapply(names(results_list), function(de_entry_name) {
    # Extract individual entry
    de_entry <- results_list[[de_entry_name]]
    de_table <- de_entry$de_table
    value_column <- de_entry$value_column
    feature_column <- de_entry$feature_column

    # Handle rownames as feature column
    if (feature_column == "rownames") {
      de_table[[feature_column]] <- rownames(de_table)
    }

    de_table[[feature_column]] <- remove_kegg_prefix(de_table[[feature_column]])

    data.frame(
      KEGG = de_table[[feature_column]], plot_value = de_table[[value_column]],
      source = rep(de_entry_name, nrow(de_table)), stringsAsFactors = FALSE
    )
  })

  return(do.call(rbind, results))
}


#' Add color palettes
#' @param nodes_df Data frame of nodes with 'value' and 'source' columns.
#' @param palettes A vector of color palette names from RColorBrewer.
#' @return nodes_df with colored nodes based on their values.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
#' @noRd
add_colors_to_nodes <- function(nodes_df, palettes = c("RdBu")) {
  # Get unique sources
  sources <- unique(na.omit(nodes_df$source))
  valid_nodes <- nodes_df[!is.na(nodes_df$source), , drop = FALSE]

  # Apply a color palette to each source
  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- rev(brewer.pal(n = 11, name = palette)) # reverse the palette BuRe
    palette_ramp <- colorRampPalette(palette_colors)
    nodes_to_color <- valid_nodes[valid_nodes$source == sources[source_index], ,
      drop = FALSE
    ]

    if (nrow(nodes_to_color) > 1) {
      range_val <- max(abs(as.numeric(nodes_to_color$plot_value)))
    } else if (nrow(nodes_to_color) == 1) {
      range_val <- abs(as.numeric(nodes_to_color$plot_value[[1]]))
    } else {
      next
    }
    # For ggplot:
    # https://stackoverflow.com/questions/79132520/symmetric-colorbar-for-values-but-print-colorbar-for-actual-observed-values

    # Cut values using real numeric range (from -range_val to +range_val)
    # Generate breaks
    breaks_seq <- seq(-range_val, range_val, length.out = 101)

    # Use the cut function to assign colors based on the breaks (cut
    # retrieves the index of the corresponding bin) cut:
    # https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/cut
    # may be useful to add general info in the dataframe:
    # https://stackoverflow.com/questions/42217741/how-do-i-add-an-attribute-to-an-r-data-frame-while-im-making-it-with-a-function
    nodes_to_color$color <- palette_ramp(100)[as.numeric(cut(as.numeric(nodes_to_color$plot_value),
      breaks = breaks_seq, include.lowest = TRUE
    ))]

    # Update main data frame
    valid_nodes$color[valid_nodes$source == sources[source_index]] <- nodes_to_color$color
  }
  nodes_df$color[!is.na(nodes_df$source)] <- valid_nodes$color

  return(nodes_df)
}


#' Add gene names to gene nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @return Updated nodes data frame with gene names added to gene nodes.
#' @noRd
add_gene_names <- function(nodes_df) {
  # find rows that are genes (logical index)
  idx <- which(!is.na(nodes_df$type) & nodes_df$type == "gene")
  if (length(idx) == 0) {
    return(nodes_df)
  }

  # Extraction of graphic_name, handle NA
  graphics_name <- as.character(nodes_df$graphics_name)
  graphics_name[is.na(graphics_name)] <- "" # Na replaced by empty
  labels <- gsub(",.*", "", graphics_name[idx]) # take first
  labels <- trimws(labels)

  nodes_df$label[idx] <- labels
  return(nodes_df)
}

#' Add compound names to compound nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @param bfc BiocFileCache object for caching KEGG compound mappings.
#' @return Updated nodes data frame with compound names added to compound nodes.
#' @importFrom BiocFileCache BiocFileCache
#' @noRd
add_compound_names <- function(nodes_df, bfc) {
  idx <- which(!is.na(nodes_df$type) & nodes_df$type == "compound")

  if (length(idx) == 0) {
    return(nodes_df)
  }

  compounds_in_graph <- as.character(nodes_df$KEGG)
  compounds_in_graph[is.na(compounds_in_graph)] <- ""
  compounds_in_graph <- compounds_in_graph[idx]

  compounds <- get_kegg_db(bfc, "compound") # expect named vector mapping KEGG id -> name
  glycan <- get_kegg_db(bfc, "glycan") # expect named vector mapping KEGG id -> name

  # safe lookup: if not found, use original id or empty string
  labels <- vapply(compounds_in_graph, function(id) {
    val <- NA_character_
    if (grepl("^C", id)) {
      tmp <- compounds[compounds[[1]] == id, 2]
      val <- if (length(tmp) > 0) tmp[1] else NA_character_
    } else if (grepl("^G", id)) {
      tmp <- glycan[glycan[1] == id, 2]
      val <- if (length(tmp) > 0) tmp[1] else NA_character_
    }
    if (is.na(val)) {
      return(id)
    }

    val <- gsub(";.*", "", val) # take first name before ';'
    return(as.character(val))
  }, character(1))

  nodes_df$label[idx] <- labels
  return(nodes_df)
}
