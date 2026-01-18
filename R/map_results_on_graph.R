
#' Add results from combined results data frame to nodes data frame
#' @param nodes_df Data frame of nodes with a column 'ids_for_mapping'
#' @param results_combined Data frame with combined results containing columns:
#' ids_for_mapping, de_value, de_source
#' @return Updated nodes data frame with added columns: de_value, color, de_source, text
#' @noRd
add_results_nodes <- function(nodes_df, results_combined) {

  # If results is empty then return the original df
  if (is.null(results_combined)) {
    return(nodes_df)
  }

  warn <- FALSE
  # nodes_df has this structure:
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

  missing_nodes <- setdiff(required_nodes, names(nodes_df))
  missing_results <- setdiff(required_results, names(results_combined))

  if (length(missing_nodes) > 0) {
    stop("Missing columns in nodes_df: ", paste(missing_nodes, collapse = ", "))
  }
  if (length(missing_results) > 0) {
    stop("Missing columns in results_combined: ", paste(missing_results, collapse = ", "))
  }

  # Subset nodes to those with valid ids_for_mapping
  ids_nodes_to_check <- which(
    !is.na(nodes_df$ids_for_mapping) & nodes_df$ids_for_mapping != ""
  )
  nodes_to_check <- nodes_df[ids_nodes_to_check, , drop = FALSE]

  # If no nodes to check, return original
  if (nrow(nodes_to_check) == 0) {
    warning("No nodes with valid 'ids_for_mapping' found.")
    return(nodes_df)
  }

  # Explode ids_for_mapping per node
  mapping <- do.call(
    rbind,
    lapply(seq_len(nrow(nodes_to_check)), function(i) {
      data.frame(
        name = nodes_to_check$name[i],
        matched_id = strsplit(
          nodes_to_check$ids_for_mapping[i], ";", fixed = TRUE
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

  if (nrow(mapping) == 0) {
    return(nodes_df)
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
    paste0,
    collapse = ""
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
  }

  # Update original dataframe
  nodes_df$de_value[ids_nodes_to_check] <- nodes_to_check$de_value
  nodes_df$de_source[ids_nodes_to_check] <- nodes_to_check$de_source
  nodes_df$text[ids_nodes_to_check] <- nodes_to_check$text

  return(nodes_df)
}

#' Combine multiple differential expression results into a single data frame
#' @param results_list A named list where each element is a differential
#' expression result
#' @return A combined data frame with columns:
#' ids_for_mapping, de_value, de_source
#' @noRd
combine_results_in_dataframe <- function(results_list) {

  if (is.null(results_list) || length(results_list) == 0) {
    return(NULL)
  }

  results <- lapply(names(results_list), function(de_entry_name) {

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

  do.call(rbind, results)
}


#' Add color palettes
#' @param nodes_df Data frame of nodes with 'de_value' and 'de_source' columns.
#' @param palettes A vector of color palette names from RColorBrewer.
#' @importFrom RColorBrewer brewer.pal
#' @return nodes_df with colored nodes based on their values.
#' @noRd
add_colors_to_nodes <- function(nodes_df, palettes = c("RdBu")) {
  sources <- unique(na.omit(nodes_df$de_source))
  valid_nodes <- nodes_df[!is.na(nodes_df$de_source), , drop = FALSE]

  for (source_index in seq_along(sources)) {
    palette <- palettes[[((source_index - 1) %% length(palettes)) + 1]]
    palette_colors <- rev(RColorBrewer::brewer.pal(n = 11, name = palette))
    palette_ramp <- colorRampPalette(palette_colors)

    nodes_to_color <- valid_nodes[
      valid_nodes$de_source == sources[source_index], ,
      drop = FALSE
    ]

    if (nrow(nodes_to_color) > 1) {
      range_val <- max(abs(as.numeric(nodes_to_color$de_value)), na.rm = TRUE)
    } else if (nrow(nodes_to_color) == 1) {
      range_val <- abs(as.numeric(nodes_to_color$de_value[[1]]))
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

    nodes_to_color$vertex.color <- palette_ramp(100)[
      as.numeric(cut(
        as.numeric(nodes_to_color$de_value),
        breaks = breaks_seq,
        include.lowest = TRUE
      ))
    ]

    valid_nodes$vertex.color[
      valid_nodes$de_source == sources[source_index]
    ] <- nodes_to_color$vertex.color
  }

  nodes_df$vertex.color[!is.na(nodes_df$de_source)] <- valid_nodes$vertex.color
  return(nodes_df)
}
