#' Add results from combined results data frame to nodes data frame
#' @param vertices_df Data frame of nodes with a column 'ids_for_mapping'
#' @param results_combined Data frame with combined results containing columns:
#' ids_for_mapping, de_value, de_source
#' @param verbose Logical, if TRUE, prints messages about the process.
#' @return Updated nodes data frame with added columns: de_value, color, de_source, text
#' @noRd
add_results_nodes <- function(vertices_df, results_combined, verbose = FALSE) {
  # If results is empty then return the original df
  if (is.null(results_combined)) {
    if (verbose) {
      message("results_combined is NULL; returning original vertices_df.")
    }
    return(vertices_df)
  }

  warn <- FALSE
  # vertices_df has this structure:
  # ids_for_mapping  | KEGG (id is the node id)
  # 1      1111
  # 2      2222
  # 2      3333;2222
  # 3    K00001
  # 4  tst00002
  # 5 undefined
  # 6      1111

  # results_combined has this structure:
  # ids_for_mapping     | de_value | source
  # 1111        2.5        de1
  # 2222       -1.2        de2
  # 3333        0.5        de1
  # Ensure required columns exist
  required_nodes <- c(
    "name", "ids_for_mapping", "de_value", "de_source", "text", "type"
  )
  required_results <- c(
    "ids_for_mapping", "de_value", "de_source"
  )

  missing_nodes <- setdiff(required_nodes, names(vertices_df))
  missing_results <- setdiff(required_results, names(results_combined))

  if (length(missing_nodes) > 0) {
    stop("Missing columns in vertices_df: ", paste(missing_nodes, collapse = ", "))
  }
  if (length(missing_results) > 0) {
    stop("Missing columns in results_combined: ", paste(missing_results, collapse = ", "))
  }

  # Subset nodes to those with valid ids_for_mapping
  ids_nodes_to_check <- which(
    !is.na(vertices_df$ids_for_mapping) & vertices_df$ids_for_mapping != ""
  )
  nodes_to_check <- vertices_df[ids_nodes_to_check, , drop = FALSE]

  # If no nodes to check, return original
  if (nrow(nodes_to_check) == 0) {
    warning("No nodes with valid 'ids_for_mapping' found.")
    return(vertices_df)
  }

  if (verbose) {
    message("Found ", nrow(nodes_to_check), " nodes with valid 'ids_for_mapping'.")
  }

  # Explode ids_for_mapping per node
  mapping <- do.call(
    rbind,
    lapply(seq_len(nrow(nodes_to_check)), function(i) {
      data.frame(
        name = nodes_to_check$name[i],
        matched_id = strsplit(
          nodes_to_check$ids_for_mapping[i], ";",
          fixed = TRUE
        )[[1]],
        stringsAsFactors = FALSE
      )
    })
  )

  # Join results onto exploded mapping
  mapping <- merge(
    mapping,
    results_combined,
    by.x = "matched_id",
    by.y = "ids_for_mapping",
    all.x = FALSE,
    sort = FALSE
  )

  # Messages
  if (nrow(mapping) == 0) {
    if (verbose) {
      message("No matches found between nodes and results.")
    }
    return(vertices_df)
  } else {
    if (verbose) {
      message("Found ", nrow(mapping), " matches between nodes and results.")
    }
  }

  # Detect multiple matches per node
  match_counts <- table(mapping$name)
  warn_nodes <- names(match_counts[match_counts > 1])
  warn <- length(warn_nodes) > 0

  # Assign first match only
  first_hits <- mapping[!duplicated(mapping$name), ]
  idx <- match(first_hits$name, nodes_to_check$name)

  nodes_to_check$de_value[idx] <- ifelse(
    is.na(nodes_to_check$de_value[idx]),
    first_hits$de_value,
    nodes_to_check$de_value[idx]
  )

  nodes_to_check$de_source[idx] <- ifelse(
    is.na(nodes_to_check$de_source[idx]),
    first_hits$de_source,
    nodes_to_check$de_source[idx]
  )

  # Append text for all matches
  sep <- ","
  mapping$text_append <- paste0(
    "Source: ", mapping$de_source,
    sep, "Value: ", mapping$de_value,
    sep, "Id: ", mapping$matched_id,
    ";"
  )

  text_by_node <- tapply(
    mapping$text_append,
    mapping$name,
    paste0
  )

  text_idx <- match(names(text_by_node), nodes_to_check$name)
  nodes_to_check$text[text_idx] <- paste0(
    nodes_to_check$text[text_idx],
    text_by_node
  )

  # Warn if necessary
  if (warn) {
    warning(
      "Some nodes had multiple matching IDs; ",
      "only the first match was used for de_value/de_source."
    )
    if (verbose) {
      message("Nodes with multiple matches: ", paste(warn_nodes, collapse = ", "))
    }
  }

  # Update original dataframe
  vertices_df$de_value[ids_nodes_to_check] <- nodes_to_check$de_value
  vertices_df$de_source[ids_nodes_to_check] <- nodes_to_check$de_source
  vertices_df$text[ids_nodes_to_check] <- nodes_to_check$text
  vertices_df$de_text[ids_nodes_to_check] <- format(round(as.numeric(nodes_to_check$de_value), 3), nsmall = 3)

  if (verbose) {
    message("Successfully added results to ", length(idx), " nodes.")
  }

  return(vertices_df)
}

#' Combine multiple differential expression results into a single data frame
#' @param results_list A named list where each element is a differential
#' expression result
#' @param verbose Logical, if TRUE, prints messages about the process.
#' @return A combined data frame with columns:
#' ids_for_mapping, de_value, de_source
#' @noRd
combine_results_in_dataframe <- function(results_list, verbose = FALSE) {
  if (is.null(results_list) || length(results_list) == 0) {
    if (verbose) {
      message("No results to combine; returning NULL.")
    }
    return(NULL)
  }

  results <- lapply(names(results_list), function(de_entry_name) {
    if (verbose) {
      message("Processing entry: ", de_entry_name)
    }

    de_entry <- results_list[[de_entry_name]]
    de_table <- de_entry$de_table
    value_column <- de_entry$value_column
    feature_column <- de_entry$feature_column

    if (feature_column == "rownames") {
      de_table[[feature_column]] <- rownames(de_table)
    }

    ids <- remove_kegg_prefix(de_table[[feature_column]])

    data.frame(
      ids_for_mapping = ids,
      de_value = de_table[[value_column]],
      de_source = rep(de_entry_name, nrow(de_table)),
      stringsAsFactors = FALSE
    )
  })

  combined_results <- do.call(rbind, results)

  if (verbose) {
    message("Combined results data frame created with ", nrow(combined_results), " rows.")
  }

  return(combined_results)
}

#' Add color palettes
#' @param vertices_df Data frame of nodes with 'de_value' and 'de_source' columns.
#' @param palettes A vector of color palette names from RColorBrewer.
#' @param verbose Logical indicating whether to print verbose messages.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames na.omit
#' @return vertices_df with colored nodes based on their values.
#' @noRd
add_colors_to_nodes <- function(vertices_df, palette = NULL, palette_limit = NULL, palettes_limits_list = NULL, palettes_list = NULL, verbose = FALSE) {
  sources <- unique(na.omit(vertices_df$de_source))
  valid_nodes <- vertices_df[!is.na(vertices_df$de_source), , drop = FALSE]
  legend <- list()

  # Checking the palette provided
  if (!is.null(palettes_list)) {
    if (!all(sources %in% names(palettes_list))) {
      warning(
        "Some sources do not have specified palettes defaulting to 'RdBu'.",
      )
      palettes_list <- NULL
    }
  }

  palettes_list <- if (!is.null(palettes_list)) {
    palettes_list
  } else {
    setNames(
      rep("RdBu", length(sources)),
      sources
    )
  }
  
  # Checking the palette limits provided
  if (!is.null(palettes_limits_list)) {
    if (!all(sources %in% names(palettes_limits_list))) {
      warning(
        "Some sources do not have specified palette limits. defaulting to max value",
      )
      palettes_limits_list <- NULL
    }
  }
  palettes_limits_list <- if (!is.null(palettes_limits_list)) {
    palettes_limits_list
  } else {
    setNames(
      rep(palette_limit, length(sources)),
      sources
    )
  }

  for (source_index in seq_along(sources)) {
    source_name <- sources[[source_index]]
    
    palette_limit <- palettes_limits_list[[source_name]]
    palette <- palettes_list[[source_name]]

    if (verbose) {
      message("Processing source: ", source_name)
    }
    if (is.null(palette)) {
      warning("Palette for source '", source_name, "' is NULL. Defaulting to 'RdBu'.")
      palette <- "RdBu"
    }
    if (length(palette) == 1 && !palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      warning(
        "Palette '", palette, "' is not a valid RColorBrewer palette. Defaulting to 'RdBu'."
      )
      palette <- "RdBu"
      palette_colors <- rev(RColorBrewer::brewer.pal(n = 11, name = palette))
    } else if (length(palette) == 1) {
      palette_colors <- rev(RColorBrewer::brewer.pal(n = 11, name = palette))
    } else {
      palette_colors <- palette
    }

    palette_ramp <- colorRampPalette(palette_colors)

    nodes_to_color <- valid_nodes[
      valid_nodes$de_source == source_name, ,
      drop = FALSE
    ]

    # This makes the de values of the plot limited to the range specified in palette_limits, if provided. This is useful to avoid outliers dominating the color scale. If palette_limits is NULL or does not contain the source_name, no limits will be applied.
    if (!is.null(palette_limit)) {
      range_val <- palette_limit
      if (verbose) {
        message(
          "Applying palette limits for source '", source_name,
          "': [-", range_val, ", ", range_val, "]"
        )
      }
      if (max(abs(nodes_to_color$de_value), na.rm = TRUE) > range_val) {
        warning(
          "Some de_value for source '", source_name,
          "' exceed the specified palette limits. Values will be capped at [", -range_val, ", ", range_val, "]."
        )
      }
      nodes_to_color$de_value <- pmax(pmin(nodes_to_color$de_value, range_val), -range_val)
    } else if (nrow(nodes_to_color) > 1) {
      if (verbose) {
        message(
          "Calculating color range for source '", source_name,
          "' based on de_value range: [-", round(min(nodes_to_color$de_value, na.rm = TRUE), 3),
          ", ", round(max(nodes_to_color$de_value, na.rm = TRUE), 3), "]"
        )
      }
      range_val <- max(abs(as.numeric(nodes_to_color$de_value)), na.rm = TRUE)
    } else if (nrow(nodes_to_color) == 1) {
      if (verbose) {
        message(
          "Only one node with de_value for source '", source_name,
          "'. Using absolute value of that node for color range: ", round(as.numeric(nodes_to_color$de_value[[1]]), 3)
        )
      }
      range_val <- abs(as.numeric(nodes_to_color$de_value[[1]]))
    } else {
      next
    }
    if (!is.finite(range_val)){
      warning(
        "Range value for source '", source_name,
        "' is not defined. Skipping color assignment for this source."
      )
      next
    } 

    legend[[source_name]] <- create_legend_single(
      range_val = range_val,
      palette_ramp = palette_ramp,
      title = source_name
    )
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

    nodes_to_color$vertex.color <- palette_ramp(100)  [
      as.numeric(cut(
        as.numeric(nodes_to_color$de_value),
        breaks = breaks_seq,
        include.lowest = TRUE
      ))
    ]

    valid_nodes$vertex.color[
      valid_nodes$de_source == source_name
    ] <- nodes_to_color$vertex.color

    if (verbose) {
      message(
        "Assigned colors for source '", source_name,
        "' using palette '", palette, "'."
      )
    }
  }

  # Create legend plot
  legend_plot <- cowplot::plot_grid(
    plotlist = legend,
    ncol = length(legend)
  )

  vertices_df$vertex.color[!is.na(vertices_df$de_source)] <- valid_nodes$vertex.color
  return_list <- list(
    vertices_df = vertices_df,
    legend_plot = legend_plot
  )
  return(return_list)
}

# Helper function to create a legend for a single source
create_legend_single <- function(range_val, palette_ramp, title = "Legend", n_element = 7) {
  # Reverse palette from RColorBrewer

  # Compute breaks based on the range of the values
  breaks_seq <- round(seq(-range_val, range_val, length.out = n_element), 1)

  # Assign colors to breaks
  legend_df <- data.frame(
    value = breaks_seq,
    color = palette_ramp(n_element)
  )

  # Create the legend plot
  p <- ggplot(legend_df) +
    geom_point(
      aes(x = 1, y = seq_along(value), fill = value),
      shape = 21,
      size = 5,
      color = "black"
    ) +
    scale_fill_gradientn(
      colours = legend_df$color,
      breaks = legend_df$value,
      name = title
    ) +
    theme_void()

  # Extract legend grob
  get_legend <- function(plot) {
    g <- ggplot_gtable(ggplot_build(plot))
    g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  }

  get_legend(p)
}
