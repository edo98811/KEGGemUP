
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
  
  nodes_to_check <- nodes_df[nodes_df$type != "line_point", , drop = FALSE]

  # Explode KEGG IDs per node 
  mapping <- do.call(
    rbind,
    lapply(seq_len(nrow(nodes_to_check)), function(i) {
      data.frame(
        id = nodes_to_check$id[i],
        KEGG = strsplit(nodes_to_check$KEGG[i], ";", fixed = TRUE)[[1]],
        stringsAsFactors = FALSE
      )
    })
  )

  # Join results onto exploded mapping
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

  # Detect multiple matches per node
  match_counts <- table(mapping$id)
  warn_nodes <- names(match_counts[match_counts > 1])
  warn <- length(warn_nodes) > 0

  # Assign first match only 
  first_hits <- mapping[!duplicated(mapping$id), ]
  idx <- match(first_hits$id, nodes_to_check$id)

  nodes_to_check$plot_value[idx] <- ifelse(
    is.na(nodes_to_check$plot_value[idx]),
    first_hits$plot_value,
    nodes_to_check$plot_value[idx]
  )

  nodes_to_check$source[idx] <- ifelse(
    is.na(nodes_to_check$source[idx]),
    first_hits$source,
    nodes_to_check$source[idx]
  )

  # Append text for all matches
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

  text_idx <- match(names(text_by_node), nodes_to_check$id)
  nodes_to_check$text[text_idx] <- paste0(
    nodes_to_check$text[text_idx],
    text_by_node
  )

  # Warn if necessary 
  if (warn) {
    warning(
      "Some nodes had multiple matching KEGG IDs; ",
      "only the first match was used for plot_value/source."
    )
  }

  nodes_df$plot_value[nodes_df$type != "line_point"] <- nodes_to_check$plot_value
  nodes_df$source[nodes_df$type != "line_point"] <- nodes_to_check$source
  nodes_df$text[nodes_df$type != "line_point"] <- nodes_to_check$text 
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
  # valid_nodes <- valid_nodes[!valid_nodes$type == "line_point", , drop = FALSE]

  # Apply a color palette to each source
  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- rev(brewer.pal(n = 11, name = palette)) # reverse the palette BuRe
    palette_ramp <- colorRampPalette(palette_colors)
    nodes_to_color <- valid_nodes[valid_nodes$source == sources[source_index], ,
      drop = FALSE
    ]

    if (nrow(nodes_to_color) > 1) {
      range_val <- max(abs(as.numeric(nodes_to_color$plot_value)), na.rm = TRUE)
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

