#' Transform a ggkegg graph to igraph or visNetwork
#'
#' @param path_id KEGG pathway ID (e.g., "hsa:04110" or "04110").
#' @param organism KEGG organism code (e.g., "hsa" for human, "mmu" for mouse).
#' @param de_results Named list of differential expression results. Each entry should be a list with elements: de_table (data.frame), value_column (character), feature_column (character), threshold (numeric).
#' @param return_type Output type: "igraph" or "visNetwork".
#' @param scaling_factor Numeric factor to scale node sizes.
#' @return An igraph or visNetwork object representing the pathway.
#'
#' @importFrom igraph graph_from_data_frame graph_attr make_empty_graph add_vertices delete_edges E V
#'
#' @export
kegg_to_graph <- function(path_id,
                          organism = "mmu",
                          de_results = NULL,
                          return_type = "igraph",
                          scaling_factor = 1.5) {
  # --- 0. Validate inputs ---
  if (!is_valid_pathway(path_id)) stop("Invalid KEGG pathway ID format.")
  if (!is.character(organism) || length(organism) != 1) stop("organism must be a single character string")
  organism <- to_organism_kegg(organism)

  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)

  # --- 1. Download KGML ---
  kgml_file <- download_kgml(path_id, bfc_kegg)
  if (is.null(kgml_file)) {
    warning("Failed to download KGML file for pathway ID: ", path_id)
    return(NULL)
  }

  # --- 2. Parse nodes and edges ---
  nodes_df <- parse_kgml_entries(kgml_file)
  edges_df <- parse_kgml_relations(kgml_file)
  nodes_df <- add_columns_nodes_df(nodes_df)

  # Extract KEGG IDs
  nodes_df$KEGG <- vapply(nodes_df$kegg_name, remove_kegg_prefix_str, character(1))

  # --- 3. Style nodes and edges ---
  nodes_df <- style_nodes(nodes_df)
  nodes_df <- add_gene_names(nodes_df)
  nodes_df <- add_compound_names(nodes_df, bfc_map)
  nodes_df <- scale_dimensions(nodes_df, factor = scaling_factor)
  nodes_df <- add_tooltip(nodes_df)
  nodes_df <- add_group(nodes_df)
  nodes_df <- nodes_df[order(nodes_df$label), ]

  if (nrow(edges_df)) edges_df <- style_edges(edges_df)

  # --- 5. Build pathway name ---
  pathway_name <- paste0("(", path_id, ") ", get_pathway_name(path_id))

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
#' @param g An igraph object representing the pathway.
#' @param de_results Named list of differential expression results.
#' @param return_type Output type: "igraph" or "visNetwork".
#' @return An igraph or visNetwork object with mapped results.
#' @importFrom visNetwork visIgraph visPhysics visLegend visOptions
#' @importFrom igraph as_data_frame graph_from_data_frame graph_attr permute V E
#' @export
map_results_to_nodes <- function(g, de_results, return_type = "visNetwork") {
  message("Mapping differential expression results to nodes...")

  # Validate each entry in de_results
  if (!is.null(de_results)) {
    # If input is a data.frame, convert to default named list
    if (inherits(de_results, "data.frame")) {
      warning("Using defaults. For personalisation use a named list of de results.")
      de_results <- list(de_input = list(
        de_table = de_results,
        value_column = "log2FoldChange",
        feature_column = "KEGG_ids"
      ))
    }

    # Check that de_results is a named list
    if (!is.list(de_results) || is.null(names(de_results)) || any(names(de_results) == "")) {
      warning("de_results must be a named list or NULL. Ignoring de_results.")
      de_results <- NULL
    }

    # Keep only valid entries
    de_results <- de_results[vapply(names(de_results), function(name) check_de_entry(de_results[[name]], name), logical(1))]
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
    nodes_df <- add_colors_to_nodes(nodes_df)

    # --- 3. Color and style nodes and edges ---
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
make_vis_graph <- function(nodes_df, edges_df, pathway_name) {
  # Different handling if no edges
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = nodes_df, main = pathway_name) # if graph has no edges
  } else {
    v <- visNetwork::visNetwork(nodes = nodes_df, edges = edges_df, main = pathway_name) # if graph has no edges
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)

  v <- visNetwork::visOptions(
    v,
    highlightNearest = list(enabled = TRUE, degree = 2, hover = TRUE),
    selectedBy = "group",
    nodesIdSelection = TRUE
  )

  v <- visNetwork::visInteraction(v, dragNodes = TRUE)
  v
  return(v)
}

#' Create an igraph graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return An igraph object representing the graph.
make_igraph_graph <- function(nodes_df, edges_df, pathway_name) {
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    fake_edges <- data.frame(from = nodes_df$name[1], to = nodes_df$name[1])
    g <- igraph::graph_from_data_frame(fake_edges, directed = FALSE, vertices = nodes_df)
    g <- igraph::delete_edges(g, igraph::E(g))
    # remove the fake edge
    # g <- make_empty_graph(n = 0, directed = FALSE)
    # g <- add_vertices(g, nrow(nodes_df), attr = as.list(nodes_df))
  } else {
    g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
  }

  g <- igraph::permute(g, order(igraph::V(g)$label))
  igraph::graph_attr(g, "title") <- pathway_name
  return(g)
}

#' Add group information to nodes based on 'undefined' groups.
#' @param nodes_df Data frame of nodes with columns: id, kegg_name, components.
#' @return nodes_df with added 'group' column for visNetwork grouping.
add_group <- function(nodes_df) {
  undefined_idx <- which(!is.na(nodes_df$kegg_name) & nodes_df$kegg_name == "undefined")

  if (length(undefined_idx) == 0) {
    return(nodes_df)
  }

  undefined_nodes <- nodes_df[undefined_idx, , drop = FALSE]

  for (i in seq_len(nrow(undefined_nodes))) {
    if (is.na(undefined_nodes$components[i]) || undefined_nodes$components[i] == "") next # If group is NA
    ids <- strsplit(undefined_nodes$components[i], ";", fixed = TRUE)[[1]]
    ids <- append(ids, undefined_nodes$id[i])

    group_label <- paste0("group_", undefined_nodes$id[i])
    nodes_df[nodes_df$id %in% ids, "group"] <- group_label
  }

  return(nodes_df)
}

#' Scale node dimensions for better visualization.
#' @param nodes_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return nodes_df with scaled x and y coordinates.
scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- -as.numeric(nodes_df$y) * factor # Invert y-axis

  return(nodes_df)
}

#' Add tooltips to nodes for visNetwork visualization.
#' @param nodes_df Data frame of nodes with columns: KEGG, label, source, value.
#' @return nodes_df with added 'title' column for tooltips.
#' @details The tooltip includes a button to the specific KEGG entry page. If multiple KEGG IDs are present, they are concatenated with '+' in the URL. It also adds information about the node name, source of differential expression data, and value.
add_tooltip <- function(nodes_df) {
  base_link <- "https://www.kegg.jp/entry/"
  #  base_link, gsub(" ", "+", nodes_df$kegg_name)

  button_html <- ifelse(
    is.na(nodes_df$kegg_name) | nodes_df$kegg_name == "",
    "",
    # paste0(
    #   '<a href="', nodes_df$link,
    #   '" target="_blank" rel="noopener noreferrer" style="text-decoration:none;">',
    #   '<button type="button">KEGG entry</button></a><br>'
    # )
    paste0(
      '<div style="text-align:center; margin-top:5px;">',
      '<a href="', ifelse(is.na(nodes_df$link), nodes_df$link, NA_character_),
      '" target="_blank" rel="noopener noreferrer" style="text-decoration:none;">',
      '<button type="button" style="',
      "background-color:#4CAF50;",
      "color:white;",
      "border:none;",
      "padding:5px 12px;",
      "border-radius:5px;",
      "cursor:pointer;",
      'font-size:12px;">',
      "KEGG entry</button></a></div>"
    )
  )

  # Ensure no "NA" strings in tooltip
  safe <- function(x) ifelse(is.na(x), "", as.character(x))
  nodes_df$title <- ifelse(
    safe(nodes_df$kegg_name) == "undefined",
    paste0("Group Node Placeholder: ", nodes_df$group, "<br>"),
    paste0(
      "Name: ", safe(nodes_df$kegg_name), "<br>",
      "Source: ", safe(nodes_df$source), "<br>",
      "Value: ", safe(format(round(as.numeric(nodes_df$value), 3), nsmall = 3)), "<br>",
      button_html
    )
  )
  return(nodes_df)
}

#' Style nodes based on their type for visNetwork visualization.
#' @param nodes_df Data frame of nodes with a column 'type'.
#' @param node_size_multiplier Numeric factor to scale node sizes (default: 1.2).
#' @return nodes_df with added visual styling columns: shape, fixed, widthConstraint, heightConstraint, size.
style_nodes <- function(nodes_df, node_size_multiplier = 1.2) {
  # Base visual settings

  # nodes_df$shape <- ifelse(nodes_df$type == "compound", "dot", "rectangle")
  nodes_df$shape <- ifelse(nodes_df$type == "compound", "dot", "box")

  # Set size constraints for non-compound nodes (compute numeric vectors first)
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


#' Style edges based on their subtype for visNetwork visualization.
#' @param edges_df Data frame of edges with a column 'subtype'.
#' @return edges_df with added visual styling columns: color, dashes, arrows, label.
style_edges <- function(edges_df) {
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
    phosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+p"),
    dephosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "-p"),
    glycosylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+g"),
    ubiquitination = list(color = "black", dashes = FALSE, arrows = "to", label = "+u"),
    methylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+m"),
    others_unknown = list(color = "black", dashes = TRUE, arrows = "to", label = "?"),
    group_relation = list(color = "transparent", dashes = TRUE, arrows = "", label = "") # Invisible edges for group relations (only functional)
  )

  # https://builtin.com/data-science/and-in-r#:~:text=The%20single%20sign%20version%20%7C%20returns,first%20element%20of%20each%20vector.
  edges_df$subtype <- gsub("/", "_", edges_df$subtype)
  edges_df$subtype[is.na(edges_df$subtype) | !(edges_df$subtype %in% names(edge_style_map))] <- "others_unknown"

  edges_df$color <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$color, character(1))
  edges_df$dashes <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$dashes, logical(1))
  edges_df$arrows <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$arrows, character(1))
  edges_df$label <- vapply(edges_df$subtype, function(x) edge_style_map[[x]]$label, character(1))

  return(edges_df)
}


#' Add results from combined results data frame to nodes data frame
#' @param nodes_df Data frame of nodes with a column 'KEGG'
#' @param results_combined Data frame with combined results containing columns: KEGG, value, source
#' @return Updated nodes data frame with added columns: value, color, source, text
add_results_nodes <- function(nodes_df, results_combined) {
  # If results is empty then return the original df
  if (is.null(results_combined)) {
    return(nodes_df)
  }

  # Iterate through nodes_df and results_combined to map values
  # For each node iterate over all results_combined
  for (i in seq_len(nrow(nodes_df))) {
    for (j in seq_len(nrow(results_combined))) {
      pattern <- results_combined$KEGG[j]
      if (is.na(pattern) || pattern == "" || is.na(nodes_df$KEGG[i]) || nodes_df$KEGG[i] == "") next
      node_ids <- strsplit(nodes_df$KEGG[i], ";", fixed = TRUE)[[1]]

      if (pattern %in% node_ids) { # KEGG ids match (substring)
        # If no value assigned to the node, assign the one from results_combined
        if (is.na(nodes_df$value[i])) {
          nodes_df$value[i] <- results_combined$value[j]
          nodes_df$source[i] <- results_combined$source[j]
          # nodes_df$text[i] <- list()
        } else { # If value warn
          warning(paste0("Multiple results mapped to node ", nodes_df$kegg_name[i], ". Keeping the first occurrence to color the node."))
        }

        # This will be added in any case
        sep <- ","
        nodes_df$text[i] <- paste0(
          nodes_df$text[i],
          "Source: ", results_combined$source[j], sep,
          "Value: ", results_combined$value[j], sep,
          "Id: ", results_combined$KEGG[j], ";"
        )
      }
    }
  }

  return(nodes_df)
}

#' Combine multiple differential expression results into a single data frame
#' @param results_list A named list where each element is a differential expression result containing a data frame (de_table), value column name (value_column), and feature column name (feature_column)
#' @return A combined data frame with columns: KEGG, value, source
combine_results_in_dataframe <- function(results_list) {
  # If no results provided, return NULL
  if (is.null(results_list) || length(results_list) == 0) {
    return(NULL)
  }

  # Combine all results into a single data frame
  results <- lapply(names(results_list), function(de_entry_name) {
    de_entry <- results_list[[de_entry_name]]
    de_table <- de_entry$de_table
    value_column <- de_entry$value_column
    feature_column <- de_entry$feature_column

    de_table[[feature_column]] <- remove_kegg_prefix(de_table[[feature_column]])

    data.frame(
      KEGG = de_table[[feature_column]],
      value = de_table[[value_column]],
      source = rep(de_entry_name, nrow(de_table)),
      stringsAsFactors = FALSE
    )
  })

  return(do.call(rbind, results))
}


#' Add color palettes
#' @param nodes_df Data frame of nodes with 'value' and 'source' columns.
#' @return nodes_df with colored nodes based on their values.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.omit
#' @importFrom grDevices colorRampPalette
add_colors_to_nodes <- function(nodes_df) {
  palettes <- c(
    "RdBu",
    "PuOr",
    "RdGy",
    "PRGn",
    "PiYG",
    "BrBG"
  )

  # Get unique sources
  sources <- unique(na.omit(nodes_df$source))
  valid_nodes <- nodes_df[!is.na(nodes_df$source), , drop = FALSE]

  # Apply a color palette to each source
  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- brewer.pal(n = 11, name = palette)
    palette_ramp <- colorRampPalette(palette_colors)
    nodes_to_color <- valid_nodes[valid_nodes$source == sources[source_index], , drop = FALSE]

    if (nrow(nodes_to_color) > 1) {
      range_val <- max(abs(as.numeric(nodes_to_color$value)))
    } else if (nrow(nodes_to_color) == 1) {
      range_val <- abs(as.numeric(nodes_to_color$value[[1]]))
    } else {
      next
    }
    # For ggplot: https://stackoverflow.com/questions/79132520/symmetric-colorbar-for-values-but-print-colorbar-for-actual-observed-values

    # Cut values using real numeric range (from -range_val to +range_val)
    # Generate breaks
    breaks_seq <- seq(-range_val, range_val, length.out = 101)

    # Use the cut function to assign colors based on the breaks (cut retrieves the index of the corresponding bin)
    # cut: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/cut
    # may be useful to add general info in the dataframe: https://stackoverflow.com/questions/42217741/how-do-i-add-an-attribute-to-an-r-data-frame-while-im-making-it-with-a-function
    nodes_to_color$color <- palette_ramp(100)[
      as.numeric(cut(as.numeric(nodes_to_color$value), breaks = breaks_seq, include.lowest = TRUE))
    ]

    # Update main data frame
    valid_nodes$color[valid_nodes$source == sources[source_index]] <- nodes_to_color$color
  }
  nodes_df$color[!is.na(nodes_df$source)] <- valid_nodes$color

  return(nodes_df)
}

#' Download and cache KEGG KGML files.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param bfc BiocFileCache object for caching KEGG KGML files.
#' @return Path to the cached KGML file.
#'
#' @importFrom xml2 read_xml
#' @importFrom httr GET http_error content status_code
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew
download_kgml <- function(pathway_id, bfc) {
  # cache key / name
  rname <- paste0(pathway_id, ".xml")

  # Check cache
  qr <- bfcquery(bfc, rname, field = "rname")

  if (nrow(qr) > 0) {
    message("Using cached KEGG KGML for ", pathway_id)
    return(bfcpath(bfc, qr$rid[1]))
  }

  # If not cached, download
  message("Downloading KGML for ", pathway_id, "...")

  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")

  # Function to get KGML content
  kgml_content <- get_kgml(url)

  if (is.null(kgml_content)) {
    warning("Failed to download KGML for ", pathway_id)
    return(NULL)
  }

  # create cache entry
  dest <- bfcnew(bfc, rname = rname, ext = ".xml")

  # write content
  writeBin(kgml_content, dest)

  message("Downloaded & cached: ", pathway_id)

  return(dest)
}


#' Add gene names to gene nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @return Updated nodes data frame with gene names added to gene nodes.
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
add_compound_names <- function(nodes_df, bfc) {

  idx <- which(!is.na(nodes_df$type) & nodes_df$type == "compound")

  if (length(idx) == 0) {
    return(nodes_df)
  }

  compounds_in_graph <- as.character(nodes_df$KEGG)
  compounds_in_graph[is.na(compounds_in_graph)] <- ""
  compounds_in_graph <- compounds_in_graph[idx]

  compounds <- get_kegg_compounds(bfc) # expect named vector mapping KEGG id -> name
  glycan <- get_kegg_glycans(bfc) # expect named vector mapping KEGG id -> name

  # safe lookup: if not found, use original id or empty string
  labels <- vapply(compounds_in_graph, function(id) {
    val <- NA_character_
    if (grepl("^C", id)) {
      val <- compounds[id]
    } else if (grepl("^G", id)) {
      val <- glycan[id]
    }
    if (is.null(val) || is.na(val)) {
      return(id)
    }
    
    val <- gsub(";.*", "", val) # take first name before ';'
    return(as.character(val))
  }, character(1))

  nodes_df$label[idx] <- labels
  return(nodes_df)
}
