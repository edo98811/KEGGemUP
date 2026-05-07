# mapping values ----------------------------------------------------------

#' Add results from combined results data frame to nodes data frame
#'
#' @param vertices_df Data frame of nodes with a column 'ids_for_mapping'
#' @param results_combined Data frame with combined results containing columns:
#' ids_for_mapping, de_value, de_source
#' @param verbose Logical, if TRUE, prints messages about the process.
#'
#' @details
#'   Example structure:
#'   - vertices_df:
#'     | ids_for_mapping |
#'     |----------------|
#'     | 1111           |
#'     | 2222           |
#'     | 3333;2222      |
#'     | K00001         |
#'     | tst00002       |
#'     | 1111           |
#'   - results_combined:
#'     | ids_for_mapping | de_value | de_source |
#'     |-----------------|----------|-----------|
#'     | 1111            | 2.5      | de1       |
#'     | 2222            | -1.2     | de2       |
#'     | 3333            | 0.5      | de1       |
#'
#' @return Updated nodes data frame with added columns: de_value, color, de_source, text
#'
#' @noRd
add_results_nodes <- function(vertices_df,
                              results_combined,
                              verbose = FALSE) {
  # If results is empty then return the original df
  if (is.null(results_combined)) {
    if (verbose) {
      message("results_combined is NULL; returning original vertices_df.")
    }
    return(vertices_df)
  }

  # Checks input validity
  if (!all(c(
    !is.null(vertices_df),
    all(c("name", "ids_for_mapping") %in% colnames(vertices_df)),
    all(c("ids_for_mapping", "de_value", "de_source") %in% colnames(results_combined)),
    dataframe_columns_are_of_type(
      vertices_df,
      c(
        "ids_for_mapping" = "character"
      ),
      # Does not accept NA values in these case.
      # They should all be either mapped or ""
      na_ok = FALSE
    )
  ))) {
    stop("add_results_nodes: Invalid input")
  }

  # To warn in case of more matches per node
  warn <- FALSE

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
    message(
      "Found ", nrow(nodes_to_check),
      " nodes with valid 'ids_for_mapping'."
    )
  }

  ### MAIN LOGIC
  # Create mapping data frame by exploding nodes with multiple ids_for_mapping
  mapping <- make_mapping_df(nodes_to_check)

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

  # Update de_value and de_source for the first match
  # (if multiple matches, only the first will be used for these columns)
  nodes_to_check$de_value[idx] <- ifelse(
    is.na(nodes_to_check$de_value[idx]),
    first_hits$de_value,
    nodes_to_check$de_value[idx]
  )

  # Update de_source for the first match
  # (if multiple matches, only the first will be used for this column)
  # If de source is NA, then assign the de source, otherwise keep it empty
  nodes_to_check$de_source[idx] <- ifelse(
    is.na(nodes_to_check$de_source[idx]),
    first_hits$de_source,
    nodes_to_check$de_source[idx]
  )

  nodes_to_check$de_name[idx] <- ifelse(
    is.na(nodes_to_check$de_name[idx]),
    first_hits$de_name,
    nodes_to_check$de_name[idx]
  )

  # Update text
  # Build the 'text' field for nodes:
  # 1. For each matched ID of a node, create a small string summarizing that match
  # 2. Group all strings by node name, concatenating them together so each node’s text
  #    includes info for all its matched IDs.
  # 3. Append this combined string to the node’s existing 'text' field.
  # Result: the 'text' field contains a complete summary of all matches for that node,
  # while numeric fields like de_value/de_source still only the first match.
  sep <- ","
  mapping$de_text_append <- paste0(
    "Source: ", mapping$de_source,
    sep, "Value: ", mapping$de_value,
    sep, "Id: ", mapping$matched_id,
    ";"
  )

  text_by_node <- tapply(
    mapping$de_text_append,
    mapping$name,
    paste0
  )

  text_idx <- match(names(text_by_node), nodes_to_check$name)
  nodes_to_check$de_text[text_idx] <- paste0(
    nodes_to_check$de_text[text_idx],
    text_by_node
  )

  # Warn if necessary
  if (warn) {
    warning(
      "Some nodes had multiple matching IDs; ",
      "only the first match was used for de_value/de_source."
    )
    if (verbose) {
      message(
        "Nodes with multiple matches: ",
        paste(warn_nodes, collapse = ", ")
      )
    }
  }

  # Update original dataframe
  vertices_df$de_value[ids_nodes_to_check] <- nodes_to_check$de_value
  vertices_df$de_source[ids_nodes_to_check] <- nodes_to_check$de_source
  vertices_df$de_text[ids_nodes_to_check] <- nodes_to_check$de_text
  vertices_df$de_name[ids_nodes_to_check] <- nodes_to_check$de_name

  if (verbose) {
    message("Successfully added results to ", length(idx), " nodes.")
  }

  return(vertices_df)
}


#' Add color palettes
#'
#' @param vertices_df Data frame of nodes with 'de_value' and 'de_source' columns.
#' @param palettes A vector of color palette names from RColorBrewer.
#' @param verbose Logical indicating whether to print verbose messages.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames na.omit
#'
#' @return vertices_df with colored nodes based on their values.
#'
#' @noRd
add_colors_to_nodes <- function(vertices_df,
                                palettes_limits_list,
                                palettes_list,
                                palette,
                                palette_limit,
                                verbose = FALSE) {
  # Default values
  sources <- unique(na.omit(vertices_df$de_source))
  valid_nodes <- vertices_df[!is.na(vertices_df$de_source), , drop = FALSE]
  legend <- list()
  default_palette <- "RdBu"

  # Checking the palette provided
  palettes_limits_list <- validate_palette_limits(
    sources = sources,
    palettes_limits_list = palettes_limits_list,
    palette_limit = palette_limit,
    default_palette_limit = FALSE
  )
  palettes_list <- validate_palettes(
    sources = sources,
    palettes_list = palettes_list,
    palette = palette,
    default_palette = default_palette
  )

  for (source_index in seq_along(sources)) {
    # Extract source name and corresponding palette/limit
    source_name <- sources[[source_index]]
    current_palette_limit <- palettes_limits_list[[source_name]]
    current_palette <- palettes_list[[source_name]]

    # Get valid nodes
    nodes_to_color <- valid_nodes[
      valid_nodes$de_source == source_name, ,
      drop = FALSE
    ]

    ## from here to map to continous value (make function)
    # Get paletteRamp for this source
    palette_ramp <- tryCatch(
      {
        colorRampPalette(current_palette)
      },
      error = function(e) {
        warning(
          "Failed to create color ramp for source '", source_name,
          "': ", e$message, ". Skipping this source."
        )
        NULL
      }
    )

    # If palette ramp creation failed, skip to next source
    if (is.null(palette_ramp)) next

    # This makes the de values of the plot limited to the range
    # specified in palette_limits, if provided.
    range_val <- get_palette_range(
      nodes_to_color$de_value,
      source_name,
      current_palette_limit,
      verbose = verbose
    )
    if (is.na(range_val)) {
      warning(
        "Could not determine a valid range for source '", source_name,
        "'. Skipping color assignment for this source."
      )
      next
    }

    # Create legend for this source
    legend[[source_name]] <- create_legend_continous(
      range_val = range_val,
      palette_ramp = palette_ramp,
      title = source_name
    )

    breaks_seq <- seq(-range_val, range_val, length.out = 101)

    # Assign colors to nodes based on their de_value
    nodes_to_color$vertex.color <- palette_ramp(100)[
      as.numeric(cut(
        winsorize(as.numeric(nodes_to_color$de_value), range_val),
        breaks = breaks_seq,
        include.lowest = TRUE
      ))
    ]
    ### end mapping to continuous value

    # here I can implement a discrete palette
    # input add argument: palette_type = "discrete"
    # if (palette_type == "discrete")
    # check that de_value is categorical, if not warn and default to continuous
    # assign colors based on factor levels of de_value and the provided palette
    # (if not same number of levels and colors, warn and repeat colors)
    # Add function to create a discrete legend

    # Update the colors in the valid_nodes data frame
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

  # Create legend plot (now return only the list, the user will decide how to plot it)
  # legend_plot <- cowplot::plot_grid(
  #   plotlist = legend,
  #   ncol = length(legend)
  # )
  vertices_df$vertex.color[
    !is.na(vertices_df$de_source)
  ] <- valid_nodes$vertex.color

  return_list <- list(
    vertices_df = vertices_df,
    legend_plots = legend
  )

  return(return_list)
}


#' Combine multiple differential expression results into a single data frame
#'
#' @param results_list A named list where each element is a differential
#' expression result
#' @param verbose Logical, if TRUE, prints messages about the process.
#'
#' @return A combined data frame with columns: ids_for_mapping, de_value,
#' de_source
#'
#' @noRd
combine_results_in_dataframe <- function(results_list,
                                         verbose = FALSE) {
  # Checks
  if (is.null(results_list) || length(results_list) == 0) {
    if (verbose) {
      message("No results to combine; returning NULL.")
    }
    return(NULL)
  }

  # Validate each entry in the results list
  results <- lapply(names(results_list), function(de_entry_name) {
    if (verbose) {
      message("Processing entry: ", de_entry_name)
    }

    # Validate the structure of each de_entry
    de_entry <- results_list[[de_entry_name]]
    de_table <- de_entry$de_table
    value_column <- de_entry$value_column
    feature_column <- de_entry$feature_column

    if (feature_column == "rownames") {
      de_table[[feature_column]] <- rownames(de_table)
    }

    ids <- remove_kegg_prefix(de_table[[feature_column]])

    # Create a data frame for this entry
    data.frame(
      ids_for_mapping = ids,
      de_value = de_table[[value_column]],
      de_source = rep(de_entry_name, nrow(de_table)),
      de_name = rep(value_column, nrow(de_table))
    )
  })

  # Combine all results into a single data frame
  combined_results <- do.call(rbind, results)

  if (verbose) {
    message(
      "Combined results data frame created with ",
      nrow(combined_results), " rows."
    )
  }

  return(combined_results)
}

#' winsorize values
#'
#' @param x the values to winsorize
#' @param range_val the (single numeric) value to cut the original vector to
#'
#' @returns the winsorized values
#'
#' @noRd
winsorize <- function(x,
                      range_val) {
  if (!is.numeric(range_val) || length(range_val) != 1) {
    stop("'range_val' must be a single numeric value")
  }
  pmax(pmin(x, range_val), -range_val)
}


# palette management ------------------------------------------------------

#' Validate color palettes
#'
#' @param sources Vector of character strings, what sources to use for the names/
#' sources to use
#' @param palettes_list List, which palettes to use; alternatively, a character
#' if these should be repeated
#' @param palette Character string, specifies a single palette
#' @param default_palette Character, which one to use by default (defaults to
#' `RdBu` for red-to-blue (through white))
#'
#' @returns A named list of palette color values
#'
#' @noRd
validate_palettes <- function(sources,
                              palettes_list = list(NA_character_),
                              palette = NULL,
                              default_palette = "RdBu"
) {

  if (!is.list(palettes_list) && !is.character(palettes_list)) {
    stop("`palettes_list` must be a list or character vector")
  }

  # Convert character vector to list for consistent handling
  if (is.character(palettes_list)) {
    palettes_list <- as.list(palettes_list)
  }

  names_present <- !is.null(names(palettes_list))
  names_match   <- names_present && all(sources %in% names(palettes_list))

  palettes_all_na <- all(vapply(
    palettes_list,
    function(x) length(x) == 1 && is.na(x),
    logical(1)
  ))

  palettes_valid <- !palettes_all_na && names_match

  if (!palettes_valid && names_present && !names_match) {
    warning("Names of `palettes_list` do not match `sources`.")
  }

  if (palettes_valid) {
    palettes_to_use <- palettes_list[sources]

  } else if (!is.null(palette)) {
    palettes_to_use <- rep(list(palette), length(sources))

  } else {
    if (!palettes_all_na) {
      warning(
        "Some sources do not have specified palettes.\n",
        "Defaulting all to '", default_palette, "'."
      )
    }

    palettes_to_use <- rep(list(default_palette), length(sources))
  }

  palettes_list_colors <- setNames(
    lapply(palettes_to_use, function(p) {
      if (is.character(p) && length(p) == 1) {
        get_palette_colors(p)
      } else {
        p
      }
    }),
    sources
  )

  palettes_list_colors
}

#' Validating palette limits
#'
#' @param sources Vector of character strings, what sources to use for the names/
#' sources to use
#' @param palettes_limits_list List for the palette limits
#' @param palette_limit Single numeric value
#' @param default_palette_limit Value to cap by default the palette limit
#'
#' @returns The `palettes_limits_list` after validation
#'
#' @noRd
validate_palette_limits <- function(sources,
                                    palettes_limits_list,
                                    palette_limit,
                                    default_palette_limit) {
  if (!is.numeric(palettes_limits_list)) {
    stop("`palettes_limits_list` must be a numeric vector")
  }

  if (is.null(names(palettes_limits_list)) ||
    !all(sources %in% names(palettes_limits_list))) {
    if (!all(is.na(palettes_limits_list))) {
      warning(
        "Some sources do not have specified palette limits.\n",
        "Not using limits for any source."
      )
    }

    palettes_limits_list <- c(NA_real_)
  }

  if (is.null(palette_limit)) {
    palette_limit <- default_palette_limit
  } else if (!is.numeric(palette_limit) || length(palette_limit) != 1) {
    stop("`palette_limit` must be a single numeric value")
  }

  palettes_limits_list <- if (all(is.na(palettes_limits_list))) {
    setNames(rep(palette_limit, length(sources)), sources)
  } else {
    palettes_limits_list
  }

  palettes_limits_list
}




#' get palette colors
#'
#' Helper function to get palette colors for a source
#'
#' @param palette Single string, should be a RColorBrewer palette name
#' @param default_palette Character, which one to use by default (defaults to
#' `RdBu` for red-to-blue (through white))
#' @param verbose Logical, defines verbosity of function
#'
#' @returns The values of the palette colors
#'
#' @noRd
get_palette_colors <- function(palette,
                               default_palette = "RdBu",
                               verbose = FALSE) {
  # If palette is a single name
  if (length(palette) == 1) {
    # Check if it is a valid RColorBrewer palette
    if (!palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      warning(
        "Palette '", palette, "' is not a valid RColorBrewer palette. ",
        "Defaulting to '", default_palette, "'."
      )
      palette <- default_palette
    }

    # Get palette colors from RColorBrewer (7 colors, reversed)
    palette_colors <- rev(RColorBrewer::brewer.pal(n = 7, name = palette))
  } else {
    # If palette is already a vector of colors, use it directly
    palette_colors <- palette
  }

  return(palette_colors)
}

#' Get a palette range
#'
#' Helper function to determine color range for a source
#'
#' @param de_value_vector Vector of the continuous values to map to the nodes
#' @param source_name Character string
#' @param palette_limit Palette limits definition
#' @param verbose Logical, defines verbosity of function
#'
#' @returns The range value for the palette
#'
#' @noRd
get_palette_range <- function(de_value_vector,
                              source_name,
                              palette_limit = FALSE,
                              verbose = FALSE) {
  # If a palette limit is provided ( I use false because NULL cannot be in a list)
  if (!isFALSE(palette_limit)) {
    range_val <- palette_limit
    if (verbose) {
      message(
        "Applying palette limits for source '", source_name,
        "': [-", range_val, ", ", range_val, "]"
      )
    }
    if (max(abs(de_value_vector), na.rm = TRUE) > range_val) {
      warning(
        "Some de_value for source '", source_name,
        "' exceed the specified palette limits.
        Values will be capped at [", -range_val, ", ", range_val, "]."
      )
    }

    # If no limit, determine range from data
  } else if (length(de_value_vector) > 1) {
    range_val <- max(abs(as.numeric(de_value_vector)), na.rm = TRUE)
    if (verbose) {
      message(
        "Calculating color range for source '", source_name,
        "' based on de_value range: [",
        round(min(de_value_vector, na.rm = TRUE), 3),
        ", ",
        round(max(de_value_vector, na.rm = TRUE), 3),
        "]"
      )
    }
  } else if (length(de_value_vector) == 1) {
    range_val <- abs(as.numeric(de_value_vector[[1]]))
    if (verbose) {
      message(
        "Only one node with de_value for source '", source_name,
        "'. Using absolute value of that node for color range: ",
        round(range_val, 3)
      )
    }
  } else {
    # No nodes to process
    return(NA_real_)
  }

  # Check if range_val is finite
  if (!is.finite(range_val)) {
    warning(
      "Range value for source '", source_name,
      "' is not defined. Skipping color assignment for this source."
    )
    return(NA_real_)
  }

  return(range_val)
}



# creating a legend ------------------------------------------------

#' Helper function to create a legend for a single source
#'
#' @param range_val Numeric, the maximum absolute value for the legend range.
#' @param palette_ramp A color ramp function created by colorRampPalette.
#' @param title Character, the title for the legend.
#' @param n_element Integer, the number of elements (breaks) to show in the legend.
#'
#' @importFrom ggplot2 ggplot geom_point aes scale_fill_gradientn theme_void
#' ggplot_gtable ggplot_build
#' @importFrom rlang .data
#'
#' @return A ggplot grob object representing the legend
#'
#' @noRd
create_legend_continous <- function(range_val,
                                    palette_ramp,
                                    title = "Legend",
                                    n_element = 7) {
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
      aes(x = 1, y = seq_along(.data$value), fill = .data$value),
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
  g <- ggplot_gtable(ggplot_build(p))
  legend <- g$grobs[[
    which(vapply(g$grobs,
                 function(x) x$name,
                 FUN.VALUE = character(1)) == "guide-box")
  ]]

  return(legend)
}


