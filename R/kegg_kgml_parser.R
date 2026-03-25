# https://xml2.r-lib.org/reference/xml_find_all.html matches xpath expressions
# https://www.w3.org/TR/xpath-31/ section B.2 xpath syntax

#' Map KEGG-styled nodes to visNetwork attributes
#' @param xml XML document object representing the KGML pathway
#' @param defaults A list of default node attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return vertices_df Data frame of nodes with visNetwork-compatible styling columns
#' @noRd
parse_kgml_nodes <- function(xml, defaults, verbose = FALSE) {
  # Find all entries that are not group or line (line is an attribute in graphics)
  nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type and not(@type='group') and (not(graphics) or not(graphics/@type='line'))]"
  )

  # Process each entry
  nodes_list <- lapply(nodes, function(node) {
    graphics_nodes <- xml2::xml_find_all(node, ".//graphics")
    n_rows <- max(length(graphics_nodes), 1) # at least one row per entry

    # Pre-allocate a data.frame for this entry
    entry_vertices_df <- data.frame(lapply(defaults, rep, each = n_rows))

    # Fill static attributes from entry
    entry_vertices_df$KEGG <- xml2::xml_attr(node, "name")
    entry_vertices_df$type <- xml2::xml_attr(node, "type")
    entry_vertices_df$link <- xml2::xml_attr(node, "link")
    entry_vertices_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes
    if (length(graphics_nodes) > 0) {
      for (i in seq_along(graphics_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_nodes: More graphics nodes than pre-allocated rows.")
        }

        g <- graphics_nodes[i]
        if (length(graphics_nodes) > 1) {
          entry_vertices_df$name[i] <- paste0(xml2::xml_attr(node, "id"), "_", i) # name must be unique
        } else {
          entry_vertices_df$name[i] <- as.character(xml2::xml_attr(node, "id"))
        }
        entry_vertices_df$graphics_name[i] <- xml2::xml_attr(g, "name")
        entry_vertices_df$x[i] <- as.integer(xml2::xml_attr(g, "x"))
        entry_vertices_df$y[i] <- as.integer(xml2::xml_attr(g, "y"))
        entry_vertices_df$graphics_type[i] <- xml2::xml_attr(g, "type")
        entry_vertices_df$width[i] <- as.integer(xml2::xml_attr(g, "width"))
        entry_vertices_df$height[i] <- as.integer(xml2::xml_attr(g, "height"))
        entry_vertices_df$fgcolor[i] <- xml2::xml_attr(g, "fgcolor")
        entry_vertices_df$bgcolor[i] <- xml2::xml_attr(g, "bgcolor")
      }
    }
    entry_vertices_df
  })

  # Combine all entries into a single data frame
  df <- do.call(rbind, nodes_list)
  if (verbose) message("Parsed ", nrow(df), " nodes from KGML.")
  df
}

#' Parse group nodes from KGML XML
#' @param xml XML document object representing the KGML pathway
#' @param defaults A list of default node attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return vertices_df Data frame of group nodes
#' @noRd
parse_kgml_groups <- function(xml, defaults, verbose = FALSE) {
  # Find all group entries
  group_nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type='group']"
  )
  # Process each entry
  nodes_list <- lapply(group_nodes, function(node) {
    n_rows <- 1 # keep the group node itself

    # Pre-allocate a data.frame for this entry
    entry_vertices_df <- data.frame(lapply(defaults, rep, each = n_rows))


    # Fill static attributes from entry
    entry_vertices_df$KEGG <- xml2::xml_attr(node, "name")
    entry_vertices_df$type <- xml2::xml_attr(node, "type")
    entry_vertices_df$link <- xml2::xml_attr(node, "link")
    entry_vertices_df$graphics_type <- "group"
    entry_vertices_df$reaction <- xml2::xml_attr(node, "reaction")
    entry_vertices_df$name <- as.character(xml2::xml_attr(node, "id"))

    components_nodes <- xml2::xml_find_all(node, ".//component")
    # Fill attributes from graphics nodes
    if (length(components_nodes) > 0) {
      for (i in seq_along(components_nodes)) {
        entry_vertices_df$components <-
          paste(xml2::xml_attr(components_nodes, "id"), collapse = ";")
      }
    }

    entry_vertices_df
  })

  # Combine all entries into a single data frame
  df <- do.call(rbind, nodes_list)
  if (verbose) message("Parsed ", nrow(df), " group nodes from KGML.")
  df
}

#' Map KEGG-styled edges to visNetwork attributes
#' @param xml XML document object representing the KGML pathway
#' @param defaults A list of default edge attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return df Data frame of nodes
#' @noRd
parse_kgml_lines <- function(xml, defaults, verbose = FALSE) {
  # Find all line entries (line is an attribute in graphics)
  line_nodes <- xml2::xml_find_all(
    xml,
    ".//entry[graphics and graphics/@type='line']"
  )

  # Process each entry
  nodes_list <- lapply(line_nodes, function(node) {
    graphics_nodes <- xml2::xml_find_all(node, ".//graphics")

    # Check consistency for line entries (it should hopefully be only one graphics node)
    if (length(graphics_nodes) == 0) { # if none
      warning("parse_kgml_lines: No graphics node found for line entry.")
      n_rows <- 1
      g <- NULL
    } else {
      if (length(graphics_nodes) > 1) { # if at least one
        warning("parse_kgml_lines: More than one graphics node found for line entry; using the first one.")
      }
      g <- graphics_nodes[[1]]
      coords <- as.numeric(strsplit(xml2::xml_attr(g, "coords"), ",")[[1]])
      n_rows <- max(length(coords) / 2, 1) # at least one row per entry
    }

    # Pre-allocate a data.frame for this entry
    entry_vertices_df <- data.frame(lapply(defaults, rep, each = n_rows))

    # Fill static attributes from entry
    entry_vertices_df$line_id <- xml2::xml_attr(node, "id")
    entry_vertices_df$name <- as.character(xml2::xml_attr(node, "id"))
    entry_vertices_df$KEGG <- xml2::xml_attr(node, "name")
    entry_vertices_df$type <- xml2::xml_attr(node, "type")
    entry_vertices_df$link <- xml2::xml_attr(node, "link")
    entry_vertices_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes (if at least two nodes for line, 4 coords)
    if (n_rows > 1 && !is.null(g)) { # checking g as well, just in case
      for (i in seq_len(n_rows)) {
        entry_vertices_df$point_index[i] <- i
        entry_vertices_df$name[i] <- paste0(xml2::xml_attr(node, "id"), "_", i) # name must be unique
        entry_vertices_df$x[i] <- coords[i * 2 - 1] # x coord first of each pair (2n-1 -> odd indexes)
        entry_vertices_df$y[i] <- coords[i * 2] # y coord second of each pair (2n -> even indexes)
        entry_vertices_df$graphics_type[i] <- "line"
        entry_vertices_df$graphics_name[i] <- xml2::xml_attr(g, "name")
        entry_vertices_df$fgcolor[i] <- xml2::xml_attr(g, "fgcolor")
        entry_vertices_df$bgcolor[i] <- xml2::xml_attr(g, "bgcolor")
      }
    }

    entry_vertices_df
  })

  # Combine all entries into a single data frame
  df <- do.call(rbind, nodes_list)
  if (verbose) message("Parsed ", nrow(df), " line nodes from KGML.")
  df
}

#' Parse line edges from KGML line nodes
#' @param line_vertices_df Data frame of line nodes extracted from parse_kgml_lines
#' @param defaults A list of default edge attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return edges_df Data frame of edges created from line nodes
#' @noRd
parse_kgml_lines_edges <- function(line_vertices_df, defaults, verbose = FALSE) {
  # Handle empty input
  if (is.null(line_vertices_df) || nrow(line_vertices_df) == 0) {
    return(NULL)
  }

  # Ensure correct order (by line_id and point_index)
  line_vertices_df <- line_vertices_df[
    order(line_vertices_df$line_id, line_vertices_df$point_index),
  ]

  # Total edges = sum(points - 1) per line
  n_edges <- sum(
    vapply(
      split(line_vertices_df$point_index, line_vertices_df$line_id),
      function(x) max(length(x) - 1L, 0L),
      integer(1)
    )
  )

  # Handle case with no edges (rbind handles NULL correctly,
  # but we want a warning because it is not an expected case)
  if (n_edges == 0) {
    warning("parse_kgml_lines_edges: No edges to create from line nodes.")
    return(NULL)
  }

  # Pre-allocate data.frame for edges
  edges_df <- as.data.frame(
    lapply(defaults, function(x) rep(x, n_edges)),
  
  )

  # Fill edges by connecting consecutive points of the same line
  row <- 1L
  # Iterate over unique line IDs
  for (id in unique(line_vertices_df$line_id)) {
    # Get indices of points for this line
    idx <- which(line_vertices_df$line_id == id)
    if (length(idx) < 2) next

    # Iterate over consecutive points (l-1 edges)
    for (i in seq_len(length(idx) - 1)) {
      edges_df$from[row] <- line_vertices_df$name[idx[i]]
      edges_df$to[row] <- line_vertices_df$name[idx[i + 1]]
      edges_df$type[row] <- "line"
      edges_df$reaction_name[row] <- line_vertices_df$reaction[idx[i]]
      row <- row + 1L
    }
  }

  if (verbose) message("Parsed ", nrow(edges_df), " edges from line nodes.")
  edges_df
}

#' Parse relation edges from KGML XML
#' @param xml XML document object representing the KGML pathway
#' @param defaults A list of default edge attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return edges_df Data frame of relation edges
#' @noRd
parse_kgml_relations <- function(xml, defaults, verbose = FALSE) {
  # Find all relation entries
  rels <- xml2::xml_find_all(xml, ".//relation")

  # Process each entry
  edges_list <- lapply(rels, function(relation) {
    subtype_nodes <- xml2::xml_find_all(relation, ".//subtype") # is  0...*
    n_rows <- max(length(subtype_nodes), 1) # at least one row per entry

    # Pre-allocate a data.frame for this entry
    entry_edges_df <- data.frame(lapply(defaults, rep, each = n_rows))

    # Fill static attributes from entry
    entry_edges_df$from <- xml2::xml_attr(relation, "entry1")
    entry_edges_df$to <- xml2::xml_attr(relation, "entry2")
    entry_edges_df$type <- "relation"
    entry_edges_df$relation_type <- xml2::xml_attr(relation, "type")

    # Fill attributes from subtype nodes
    if (length(subtype_nodes) > 0) {
      for (i in seq_along(subtype_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_relations: More subtype nodes than pre-allocated rows.")
        }
        g <- subtype_nodes[i]
        entry_edges_df$relation_subtype_name[i] <- gsub("[/ ]", "_", xml2::xml_attr(g, "name")) # replace / and space with _
        entry_edges_df$relation_subtype_value[i] <- xml2::xml_attr(g, "value")
      }
    }

    entry_edges_df
  })

  # Combine all entries into a single data frame
  df <- do.call(rbind, edges_list)
  if (verbose) message("Parsed ", nrow(df), " relations from KGML.")
  df
}

#' Complete reaction edges by connecting substrates and products to reaction nodes
#' @param vertices_df Data frame of nodes with a column 'reaction' containing reaction IDs
#' @param edges_df Data frame of edges with a column 'type' indicating edge type and columns 'from' and 'to' for node IDs
#' @param defaults A list of default edge attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return edges_df Updated data frame of edges with reaction edges completed
#' @noRd
complete_kgml_reactions <- function(vertices_df, edges_df, defaults, verbose = FALSE) {
  reaction_edges <- edges_df[edges_df$type == "reaction", ]
  other_edges <- edges_df[edges_df$type != "reaction", ]

  if (is.null(reaction_edges)) {
    if (verbose) message("No reaction edges to complete.")
    return(edges_df)
  }

  # Collect all new edges in a list
  new_edges_list <- lapply(seq_len(nrow(reaction_edges)), function(i) {
    reaction_id <- reaction_edges$reaction_name[i]

    # Skip if reaction_id is missing or empty, with a warning
    if (is.na(reaction_id) || reaction_id == "") {
      if (verbose) message("Skipping edge ", i, " with missing reaction name.")
      return(NULL)
    }
    reaction_node_idx <- which(vertices_df$reaction == reaction_id & vertices_df$graphics_type != "line")

    if (length(reaction_node_idx) == 0) {
      if (verbose) message("No node found for reaction ", reaction_id, " in edge ", i)
      return(NULL)
    }

    n_rows <- length(reaction_node_idx) * 2
    entry_edges_df <- data.frame(lapply(defaults, rep, each = n_rows))

    preserve_cols <- setdiff(names(reaction_edges), c("from", "to", "reaction_type"))
    for (j in seq_along(reaction_node_idx)) {
      row_offset <- (j - 1) * 2

      # Substrate -> Reaction
      entry_edges_df[row_offset + 1, preserve_cols] <- reaction_edges[i, preserve_cols]
      entry_edges_df[row_offset + 1, "from"] <- reaction_edges$from[i]
      entry_edges_df[row_offset + 1, "to"] <- vertices_df$name[reaction_node_idx[j]]
      entry_edges_df[row_offset + 1, "reaction_type"] <- paste0("reaction_substrate_", reaction_edges$reaction_type[i])

      # Reaction -> Product
      entry_edges_df[row_offset + 2, preserve_cols] <- reaction_edges[i, preserve_cols]
      entry_edges_df[row_offset + 2, "from"] <- vertices_df$name[reaction_node_idx[j]]
      entry_edges_df[row_offset + 2, "to"] <- reaction_edges$to[i]
      entry_edges_df[row_offset + 2, "reaction_type"] <- paste0("reaction_product_", reaction_edges$reaction_type[i])
    }

    entry_edges_df
  })

  # Filter NULL elements and combine
  valid_edges <- Filter(Negate(is.null), new_edges_list)

  if (length(valid_edges) == 0) {
    if (verbose) message("No valid reaction edges to add.")
    return(other_edges)
  }

  all_new_edges <- do.call(rbind, valid_edges)
  result <- rbind(other_edges, all_new_edges)

  if (verbose) message("Completed reaction edges. Total edges: ", nrow(result))
  return(result)
}

#' Parse reaction edges from KGML XML
#' @param xml XML document object representing the KGML pathway
#' @param defaults A list of default edge attributes
#' @param verbose Logical indicating whether to print verbose messages
#' @return edges_df Data frame of reaction edges
#' @noRd
parse_kgml_reactions <- function(xml, defaults, verbose = FALSE) {

  # Find all reaction entries
  reactions <- xml2::xml_find_all(xml, ".//reaction")

  # Process each entry
  nodes_list <- lapply(reactions, function(reaction) {
    # Substrates and products
    substrates_nodes <- xml2::xml_find_all(reaction, ".//substrate") # is 1...*
    products_nodes <- xml2::xml_find_all(reaction, ".//product") # is 1...*

    n_sub <- length(substrates_nodes)
    n_prod <- length(products_nodes)

    # All possible combinations of substrates and products
    n_rows <- n_sub * n_prod

    # Pre-allocate a data.frame for this entry
    entry_edges_df <- data.frame(lapply(defaults, rep, each = n_rows))

    # Fill static attributes from entry
    entry_edges_df$reaction_id <- xml2::xml_attr(reaction, "id")
    entry_edges_df$reaction_name <- xml2::xml_attr(reaction, "name")
    entry_edges_df$reaction_type <- xml2::xml_attr(reaction, "type")
    entry_edges_df$type <- "reaction"

    # Pre-extract substrate attributes
    sub_id <- xml2::xml_attr(substrates_nodes, "id")
    sub_name <- xml2::xml_attr(substrates_nodes, "name")
    sub_alt <- vapply(
      substrates_nodes,
      function(x) xml2::xml_attr(xml2::xml_find_first(x, "./alt"), "name"),
      character(1)
    )

    # Pre-extract product attributes
    prod_id <- xml2::xml_attr(products_nodes, "id")
    prod_name <- xml2::xml_attr(products_nodes, "name")
    prod_alt <- vapply(
      products_nodes,
      function(x) xml2::xml_attr(xml2::xml_find_first(x, "./alt"), "name"),
      character(1)
    )

    # Fill all possible combinations
    row <- 1L
    for (p in seq_len(n_prod)) {
      for (s in seq_len(n_sub)) {
        entry_edges_df$from[row] <- sub_id[s]
        entry_edges_df$to[row] <- prod_id[p]

        entry_edges_df$reaction_from_name[row] <- sub_name[s]
        entry_edges_df$reaction_to_name[row] <- prod_name[p]

        entry_edges_df$reaction_alt_name_substrate[row] <- sub_alt[s]
        entry_edges_df$reaction_alt_name_product[row] <- prod_alt[p]

        row <- row + 1L
      }
    }
    entry_edges_df
  })

  # Combine all entries into a single data frame
  df <- do.call(rbind, nodes_list)
  if (verbose) message("Parsed ", nrow(df), " reactions from KGML.")
  df
}


#' Add compound names to compound nodes in the nodes data frame.
#' @param vertices_df Data frame of nodes with a column 'type' indicating node type.
#' @param bfc BiocFileCache object for caching KEGG compound mappings.
#' @param verbose Logical indicating whether to print verbose messages.
#' @return Updated nodes data frame with compound names added to compound nodes.
#' @importFrom BiocFileCache BiocFileCache
#' @noRd
add_node_labels <- function(vertices_df, bfc, verbose = FALSE) {
  # Load KEGG databases
  compounds_db <- get_kegg_db(db_name = "compound", bfc = bfc, verbose = verbose)
  glycans_db <- get_kegg_db(db_name = "glycan", bfc = bfc, verbose = verbose)
  genes_db <- get_kegg_db(db_name = "ko", bfc = bfc, verbose = verbose)
  enzymes_db <- get_kegg_db(db_name = "enzyme", bfc = bfc, verbose = verbose)

  # Convert to named lookup vectors
  compounds_lookup <- setNames(as.character(compounds_db[, 2]), compounds_db[, 1])
  glycans_lookup <- setNames(as.character(glycans_db[, 2]), glycans_db[, 1])
  genes_lookup <- setNames(as.character(genes_db[, 2]), genes_db[, 1])
  enzymes_lookup <- setNames(as.character(enzymes_db[, 2]), enzymes_db[, 1])

  # Initialize labels
  map_ids <- sub("[;].*", "", vertices_df$ids_for_mapping)
  ids <- vertices_df$KEGG
  labels <- sub("[;,].*", "", vertices_df$graphics_name)
  labels[is.na(labels)] <- map_ids[is.na(labels)]

  # Compounds
  is_c <- grepl("^cpd:C", ids)
  found_c <- compounds_lookup[map_ids[is_c]]
  na_pos <- is.na(found_c)
  found_c[na_pos] <- labels[is_c][na_pos]
  labels[is_c] <- sub("[;].*", "", found_c)
  if (verbose) message("Mapped ", sum(!na_pos & is_c), " compounds.")

  # Glycans
  is_g <- grepl("^gl:G", ids)
  found_g <- glycans_lookup[map_ids[is_g]]
  na_pos <- is.na(found_g)
  found_g[na_pos] <- labels[is_g][na_pos]
  labels[is_g] <- sub("[;].*", "", found_g)
  if (verbose) message("Mapped ", sum(!na_pos & is_g), " glycans.")

  # Genes
  is_k <- grepl("^ko:", ids)
  found_k <- genes_lookup[map_ids[is_k]]
  na_pos <- is.na(found_k)
  found_k[na_pos] <- labels[is_k][na_pos]
  labels[is_k] <- sub("[;].*", "", found_k)
  if (verbose) message("Mapped ", sum(!na_pos & is_k), " genes.")

  # Enzymes
  is_e <- grepl("^ec:", ids)
  found_e <- enzymes_lookup[map_ids[is_e]]
  na_pos <- is.na(found_e)
  found_e[na_pos] <- labels[is_e][na_pos]
  labels[is_e] <- sub("[;].*", "", found_e)
  if (verbose) message("Mapped ", sum(!na_pos & is_e), " enzymes.")

  # Assign node labels
  vertices_df$label <- labels
  vertices_df
}

#' Add reaction labels to reaction nodes in the nodes data frame.
#' @param vertices_df Data frame of nodes with a column 'reaction' containing reaction IDs.
#' @param bfc BiocFileCache object for caching KEGG reaction mappings.
#' @param verbose Logical indicating whether to print verbose messages.
#' @return Updated nodes data frame with reaction labels added to reaction nodes.
#' @noRd
add_reaction_labels <- function(vertices_df, bfc, verbose = FALSE) {
  # Load reaction database
  reactions_db <- get_kegg_db(db_name = "reaction", bfc = bfc)
  reactions_lookup <- setNames(as.character(reactions_db[, 2]), reactions_db[, 1])

  # Extract reaction IDs
  reaction_ids <- vapply(vertices_df$reaction, remove_kegg_prefix_str, character(1))
  is_r <- grepl("^R", reaction_ids)

  # Lookup labels
  found_r <- reactions_lookup[reaction_ids[is_r]]
  na_pos <- is.na(found_r)
  found_r[na_pos] <- reaction_ids[is_r][na_pos]

  # Prepare reaction labels
  reaction_labels <- reaction_ids
  reaction_labels[is_r] <- sub(";.*", "", found_r)
  reaction_labels[!is_r] <- NA_character_

  # Assign columns
  vertices_df$reaction_label <- reaction_labels
  vertices_df$reaction_link <- ifelse(
    is_r,
    paste0("https://www.kegg.jp/dbget-bin/www_bget?", reaction_ids),
    NA_character_
  )

  if (verbose) {
    message("Mapped ", sum(!na_pos), " reactions.")
  }

  vertices_df
}

#' Add group labels and coordinates to group nodes in the nodes data frame.
#' @param vertices_df Data frame of nodes with a column 'type' indicating node type
#' and a column 'components' listing component node IDs.
#' @return Updated nodes data frame with group labels and coordinates added to group nodes.
#' @noRd
add_group <- function(vertices_df, verbose = FALSE) {
  # Identify undefined nodes (group nodes)
  group_idx <- which(vertices_df$type == "group")
  if (length(group_idx) == 0) {
    return(vertices_df)
  }

  # Loop only over group nodes with non-empty components
  for (i in group_idx) {
    comps <- vertices_df$components[i]
    if (is.na(comps) || comps == "") next

    # Split components and include the group node itself
    ids <- c(strsplit(comps, ";", fixed = TRUE)[[1]], vertices_df$name[i])

    # Get indices of all nodes in this group
    node_idx <- match(ids, vertices_df$name)

    # just to be safe (should not happen)
    if (any(is.na(node_idx))) {
      warning(
        "add_group: Some component IDs not found in vertices_df: ",
        paste(ids[is.na(node_idx)], collapse = ", ")
      )
      node_idx <- node_idx[!is.na(node_idx)]
    }

    # Build group label from component labels 
    # (exclude the last one, which is the group node itself)
    comp_labels <- vertices_df$label[node_idx[-length(node_idx)]]
    group_label <- paste(comp_labels, collapse = ", ")

    # Assign group label to all nodes in this group
    vertices_df$group[node_idx] <- group_label

    # Compute average x and y coordinates of all component nodes
    avg_x <- round(mean(as.numeric(vertices_df$x[node_idx]), na.rm = TRUE))
    avg_y <- round(mean(as.numeric(vertices_df$y[node_idx]), na.rm = TRUE))

    # Assign average coordinates to the group node itself
    vertices_df$x[i] <- avg_x
    vertices_df$y[i] <- avg_y

    if (verbose) {
      message(
        "Group node '", vertices_df$name[i], "' assigned label: '",
        group_label, "' at (", avg_x, ", ", avg_y, ")"
      )
    }
  }

  vertices_df
}
