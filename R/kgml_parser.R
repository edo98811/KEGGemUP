# https://xml2.r-lib.org/reference/xml_find_all.html matches xpath expressions
# https://www.w3.org/TR/xpath-31/ section B.2 xpath syntax
parse_kgml_nodes <- function(xml, defaults) {
  # Find all entries that are not group or line
  nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type and not(@type='group' or @type='line')]"
  )

  # Process each entry
  nodes_list <- lapply(nodes, function(node) {
    graphics_nodes <- xml2::xml_find_all(node, ".//graphics")
    n_rows <- max(length(graphics_nodes), 1) # at least one row per entry

    # Pre-allocate a single-entry data frame
    single_entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    single_entry_nodes_df$name <- xml2::xml_attr(node, "id")
    single_entry_nodes_df$kegg_name <- xml2::xml_attr(node, "name")
    single_entry_nodes_df$type <- xml2::xml_attr(node, "type")
    single_entry_nodes_df$link <- xml2::xml_attr(node, "link")
    single_entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes
    if (length(graphics_nodes) > 0) {
      for (i in seq_along(graphics_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_nodes: More graphics nodes than pre-allocated rows.")
        }
        g <- graphics_nodes[i]

        single_entry_nodes_df$graphics_name[i] <- xml2::xml_attr(g, "name")
        single_entry_nodes_df$x[i] <- xml2::xml_attr(g, "x")
        single_entry_nodes_df$y[i] <- xml2::xml_attr(g, "y")
        single_entry_nodes_df$graphics_type[i] <- xml2::xml_attr(g, "type")
        single_entry_nodes_df$width[i] <- xml2::xml_attr(g, "width")
        single_entry_nodes_df$height[i] <- xml2::xml_attr(g, "height")
        single_entry_nodes_df$fgcolor[i] <- xml2::xml_attr(g, "fgcolor")
      }
    }

    single_entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

parse_kgml_groups <- function(xml, defaults) {
  # Find all entries that are not group or line
  group_nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type='group']"
  )
  # Process each entry
  nodes_list <- lapply(nodes, function(node) {
    components_nodes <- xml2::xml_find_all(node, ".//component")
    n_rows <- max(length(components_nodes), 1) # at least one row per entry

    # Pre-allocate a single-entry data frame
    single_entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    single_entry_nodes_df$name <- xml2::xml_attr(node, "id")
    single_entry_nodes_df$kegg_name <- xml2::xml_attr(node, "name")
    single_entry_nodes_df$type <- xml2::xml_attr(node, "type")
    single_entry_nodes_df$link <- xml2::xml_attr(node, "link")
    single_entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes
    if (length(components_nodes) > 0) {
      for (i in seq_along(components_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_nodes: More graphics nodes than pre-allocated rows.")
        }
        g <- components_nodes[i]
        single_entry_nodes_df$components[i] <- xml2::xml_attr(g, "id")
      }
    }

    single_entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

parse_kgml_lines <- function(xml, defaults) {
  # Find all entries that are not group or line
  group_nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type='line']"
  )

  # Process each entry
  nodes_list <- lapply(nodes, function(node) {
    coords <- as.numeric(strsplit(xml2::xml_attr(node, "line"), ",")[[1]])
    if (length(xy) < 4 || any(is.na(xy))) { # at least 2 coords (x1,y1,x2,y2)
      return(NULL)
    }
    n_rows <- max(length(coords)/2, 1) # at least one row per entry

    # Pre-allocate a single-entry data frame
    single_entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    single_entry_nodes_df$name <- xml2::xml_attr(node, "id")
    single_entry_nodes_df$kegg_name <- xml2::xml_attr(node, "name")
    single_entry_nodes_df$type <- xml2::xml_attr(node, "type")
    single_entry_nodes_df$link <- xml2::xml_attr(node, "link")
    single_entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes
    if (length(components_nodes) > 0) {
      for (i in seq_along(components_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_nodes: More graphics nodes than pre-allocated rows.")
        }
        g <- components_nodes[i]
        single_entry_nodes_df$components[i] <- xml2::xml_attr(g, "id")
      }
    }

    single_entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
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

  compounds <- get_kegg_db(db_name = "compound", bfc = bfc) # expect named vector mapping KEGG id -> name
  glycan <- get_kegg_db(db_name = "glycan", bfc = bfc) # expect named vector mapping KEGG id -> name

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
