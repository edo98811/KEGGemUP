# https://xml2.r-lib.org/reference/xml_find_all.html matches xpath expressions
# https://www.w3.org/TR/xpath-31/ section B.2 xpath syntax
parse_kgml_nodes <- function(xml, defaults) {
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
    entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    entry_nodes_df$KEGG <- xml2::xml_attr(node, "name")
    entry_nodes_df$type <- xml2::xml_attr(node, "type")
    entry_nodes_df$link <- xml2::xml_attr(node, "link")
    entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes
    if (length(graphics_nodes) > 0) {
      for (i in seq_along(graphics_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_nodes: More graphics nodes than pre-allocated rows.")
        }
        g <- graphics_nodes[i]

        if (length(graphics_nodes) > 1) {
          entry_nodes_df$name[i] <- paste0(xml2::xml_attr(node, "id"), "_", i) # name must be unique
        } else {
          entry_nodes_df$name[i] <- xml2::xml_attr(node, "id")
        }
        entry_nodes_df$graphics_name[i] <- xml2::xml_attr(g, "name")
        entry_nodes_df$x[i] <- xml2::xml_attr(g, "x")
        entry_nodes_df$y[i] <- xml2::xml_attr(g, "y")
        entry_nodes_df$graphics_type[i] <- xml2::xml_attr(g, "type")
        entry_nodes_df$width[i] <- xml2::xml_attr(g, "width")
        entry_nodes_df$height[i] <- xml2::xml_attr(g, "height")
        entry_nodes_df$fgcolor[i] <- xml2::xml_attr(g, "fgcolor")
        entry_nodes_df$bgcolor[i] <- xml2::xml_attr(g, "bgcolor")
      }
    }

    entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

parse_kgml_groups <- function(xml, defaults) {
  # Find all group entries
  group_nodes <- xml2::xml_find_all(
    xml,
    ".//entry[@type='group']"
  )
  # Process each entry
  nodes_list <- lapply(group_nodes, function(node) {
    n_rows <- 1 # keep the group node itself

    # Pre-allocate a data.frame for this entry
    entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    entry_nodes_df$kegg_entry_id <- xml2::xml_attr(node, "name")
    entry_nodes_df$KEGG <- xml2::xml_attr(node, "name")
    entry_nodes_df$type <- xml2::xml_attr(node, "type")
    entry_nodes_df$link <- xml2::xml_attr(node, "link")
    entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")
    entry_nodes_df$name <- xml2::xml_attr(node, "id") # first row is the group node itself
    components_nodes <- xml2::xml_find_all(node, ".//component")
    # Fill attributes from graphics nodes
    if (length(components_nodes) > 0) {
      for (i in seq_along(components_nodes)) {
        entry_nodes_df$components <- paste(components_nodes %>% xml2::xml_attr("id"), collapse = ";")
        # g <- components_nodes[i + 1]
        # entry_nodes_df$name[i + 1] <- paste0(xml2::xml_attr(node, "id"), "_", i) # name must be unique
        # entry_nodes_df$components[i + 1] <- xml2::xml_attr(g, "id")
      }
    }

    entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

parse_kgml_lines <- function(xml, defaults) {
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
    entry_nodes_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    entry_nodes_df$line_id <- xml2::xml_attr(node, "id")
    entry_nodes_df$name <- xml2::xml_attr(node, "id")
    entry_nodes_df$KEGG <- xml2::xml_attr(node, "name")
    entry_nodes_df$type <- xml2::xml_attr(node, "type")
    entry_nodes_df$link <- xml2::xml_attr(node, "link")
    entry_nodes_df$reaction <- xml2::xml_attr(node, "reaction")

    # Fill attributes from graphics nodes (if at least two nodes for line, 4 coords)
    if (n_rows > 1 && !is.null(g)) { # checking g as well, just in case
      for (i in seq_len(n_rows)) {
        entry_nodes_df$point_index[i] <- i
        entry_nodes_df$name[i] <- paste0(xml2::xml_attr(node, "id"), "_", i) # name must be unique
        entry_nodes_df$x[i] <- coords[i * 2 - 1] # x coord first of each pair (2n-1 -> odd indexes)
        entry_nodes_df$y[i] <- coords[i * 2] # y coord second of each pair (2n -> even indexes)
        entry_nodes_df$graphics_type[i] <- "line"
        entry_nodes_df$fgcolor[i] <- xml2::xml_attr(g, "fgcolor")
        entry_nodes_df$bgcolor[i] <- xml2::xml_attr(g, "bgcolor")
      }
    }

    entry_nodes_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

parse_kgml_lines_edges <- function(line_nodes_df, defaults) {
  # Ensure correct order (by line_id and point_index)
  line_nodes_df <- line_nodes_df[
    order(line_nodes_df$line_id, line_nodes_df$point_index),
  ]

  # Total edges = sum(points - 1) per line
  n_edges <- sum(
    vapply(
      split(line_nodes_df$point_index, line_nodes_df$line_id),
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
    stringsAsFactors = FALSE
  )

  # Fill edges by connecting consecutive points of the same line
  row <- 1L
  # Iterate over unique line IDs
  for (id in unique(line_nodes_df$line_id)) {
    # Get indices of points for this line
    idx <- which(line_nodes_df$line_id == id)
    if (length(idx) < 2) next

    # Iterate over consecutive points (l-1 edges)
    for (i in seq_len(length(idx) - 1)) {
      edges_df$from[row] <- line_nodes_df$name[idx[i]]
      edges_df$to[row] <- line_nodes_df$name[idx[i + 1]]
      edges_df$type[row] <- "line"
      row <- row + 1L
    }
  }

  edges_df
}

parse_kgml_relations <- function(xml, defaults) {
  # Find all relation entries
  rels <- xml2::xml_find_all(xml, ".//relation")

  # Process each entry
  nodes_list <- lapply(rels, function(relation) {
    subtype_nodes <- xml2::xml_find_all(relation, ".//subtype") # is  0...*
    n_rows <- max(length(subtype_nodes), 1) # at least one row per entry

    # Pre-allocate a data.frame for this entry
    entry_edges_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    entry_edges_df$from <- xml2::xml_attr(relation, "entry1")
    entry_edges_df$to <- xml2::xml_attr(relation, "entry2")
    entry_edges_df$type <- xml2::xml_attr(relation, "type")

    # Fill attributes from subtype nodes
    if (length(subtype_nodes) > 0) {
      for (i in seq_along(subtype_nodes)) {
        if (i > n_rows) {
          stop("parse_kgml_relations: More subtype nodes than pre-allocated rows.")
        }
        g <- subtype_nodes[i]
        entry_edges_df$relation_subtype_name[i] <- xml2::xml_attr(g, "name")
        entry_edges_df$relation_subtype_value[i] <- xml2::xml_attr(g, "value")
      }
    }

    entry_edges_df
  })

  # Combine all entries into a single data frame
  do.call(rbind, nodes_list)
}

#' @noRd
parse_kgml_reactions <- function(xml, defaults) {
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
    entry_edges_df <- as.data.frame(
      lapply(defaults, function(x) rep(x, n_rows)),
      stringsAsFactors = FALSE
    )

    # Fill static attributes from entry
    entry_edges_df$reaction_id <- xml2::xml_attr(reaction, "id")
    entry_edges_df$reaction_name <- xml2::xml_attr(reaction, "name")
    entry_edges_df$reaction_type <- xml2::xml_attr(reaction, "type")

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
  do.call(rbind, nodes_list)
}


#' Add compound names to compound nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @param bfc BiocFileCache object for caching KEGG compound mappings.
#' @return Updated nodes data frame with compound names added to compound nodes.
#' @importFrom BiocFileCache BiocFileCache
#' @noRd
add_node_labels <- function(nodes_df, bfc) {
  # Extract KEGG IDs
  ids <- nodes_df$KEGG_no_prefix

  # Load KEGG databases
  compounds_db <- get_kegg_db(db_name = "compound", bfc = bfc)
  glycans_db <- get_kegg_db(db_name = "glycan", bfc = bfc)
  genes_db <- get_kegg_db(db_name = "ko", bfc = bfc)

  # Convert to named lookup vectors
  compounds_lookup <- setNames(as.character(compounds_db[, 2]), compounds_db[, 1])
  glycans_lookup <- setNames(as.character(glycans_db[, 2]), glycans_db[, 1])
  genes_lookup <- setNames(as.character(genes_db[, 2]), genes_db[, 1])

  # Initialize labels
  labels <- ids

  # Compounds (Cxxxxx)
  is_c <- grepl("^C", ids)
  found_c <- compounds_lookup[ids[is_c]]
  found_c[is.na(found_c)] <- ids[is_c]
  labels[is_c] <- sub(";.*", "", found_c)

  # Glycans (Gxxxxx)
  is_g <- grepl("^G", ids)
  found_g <- glycans_lookup[ids[is_g]]
  found_g[is.na(found_g)] <- ids[is_g]
  labels[is_g] <- sub(";.*", "", found_g)

  # Genes (Kxxxxx)
  is_k <- grepl("^K", ids)
  found_k <- genes_lookup[ids[is_k]]
  found_k[is.na(found_k)] <- ids[is_k]
  labels[is_k] <- sub(";.*", "", found_k)

  # Assign node labels
  nodes_df$label <- labels
  nodes_df
}

#' Add reaction labels to reaction nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'reaction' containing reaction IDs.
#' @param bfc BiocFileCache object for caching KEGG reaction mappings.
#' @return Updated nodes data frame with reaction labels added to reaction nodes.
#' @noRd
add_reaction_labels <- function(nodes_df, bfc) {
  # Load reaction database
  reactions_db <- get_kegg_db(db_name = "reaction", bfc = bfc)
  reactions_lookup <- setNames(as.character(reactions_db[, 2]), reactions_db[, 1])

  # Extract reaction IDs
  reaction_ids <- vapply(nodes_df$reaction, remove_kegg_prefix_str, character(1))
  is_r <- grepl("^R", reaction_ids)

  # Lookup labels
  found_r <- reactions_lookup[reaction_ids[is_r]]
  found_r[is.na(found_r)] <- reaction_ids[is_r]

  reaction_labels <- reaction_ids
  reaction_labels[is_r] <- sub(";.*", "", found_r)

  # Assign columns
  nodes_df$reaction_label <- reaction_labels
  nodes_df$reaction_link <- ifelse(
    is_r,
    paste0("https://www.kegg.jp/dbget-bin/www_bget?", reaction_ids),
    NA_character_
  )

  nodes_df
}


add_group <- function(nodes_df) {
  # Identify undefined nodes (group nodes)
  group_idx <- which(nodes_df$type == "group")
  if (length(group_idx) == 0) {
    return(nodes_df)
  }

  # Loop only over group nodes with non-empty components
  for (i in group_idx) {
    comps <- nodes_df$components[i]
    if (is.na(comps) || comps == "") next

    # Split components and include the group node itself
    ids <- c(strsplit(comps, ";", fixed = TRUE)[[1]], nodes_df$name[i])

    # Get indices of all nodes in this group
    node_idx <- match(ids, nodes_df$name)

    # just to be safe (should not happen)
    if (any(is.na(node_idx))) {
      warning(
        "add_group: Some component IDs not found in nodes_df: ",
        paste(ids[is.na(node_idx)], collapse = ", ")
      )
      node_idx <- node_idx[!is.na(node_idx)]
    }

    # Build group label from component labels (exclude the last one, which is the group node itself)
    comp_labels <- nodes_df$label[node_idx[-length(node_idx)]]
    group_label <- paste(comp_labels, collapse = ";")

    # Assign group label to all nodes in this group
    nodes_df$group[node_idx] <- group_label
  }

  nodes_df
}
