#' Parse KEGG KGML files to extract relations and edges data frames.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{from}{The ID of the source entry (node) in the pathway, corresponding to the `entry1` attribute in the KGML `<relation>` tag.}
#'   \item{to}{The ID of the target entry (node) in the pathway, corresponding to the `entry2` attribute in the KGML `<relation>` tag.}
#'   \item{type}{The general type of relationship between the two entries (e.g., `'ECrel'`, `'PPrel'`, `'GErel'`, `'PCrel'`, `'maplink'`). These types describe the biological nature of the connection, such as enzyme-enzyme relation or protein-protein interaction.}
#'   \item{relation_subtype}{A more specific relation_subtype of the relation, derived from the `<relation_subtype>` child elements of the `<relation>` node (e.g., `'activation'`, `'inhibition'`, `'expression'`, `'compound'`). If no relation_subtype is defined, this will be `NA`.}
#'   \item{relation_value}{A categorical value associated with the relation_subtype, which encodes the apearence of the arrow (from the `value` attribute of `<relation_subtype>`). If not present, this will be `NA`.}
#' }
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr
#'
#' @details
#' The function reads a KEGG KGML (KEGG Markup Language) file, which encodes pathway
#' information as XML, and extracts all `<relation>` elements that describe the
#' interactions or relationships between entities in the pathway. Each `<relation>` may
#' contain one or more `<relation_subtype>` elements that provide additional details about
#' the interaction. The result is a tidy data.frame suitable for network analysis or
#' visualization, where each row represents one relation–relation_subtype pair.
#'
#' @export
parse_kgml_relations <- function(file) {
  doc <- read_xml(file)

  rels <- xml_find_all(doc, ".//relation")

  rels_list <- lapply(rels, function(rel) {
    entry1 <- xml_attr(rel, "entry1")
    entry2 <- xml_attr(rel, "entry2")
    type <- xml_attr(rel, "type")
    subnodes <- xml_find_all(rel, ".//subtype")

    if (length(subnodes) == 0) {
      data.frame(
        from = entry1,
        to = entry2,
        type = type,
        relation_subtype = NA_character_,
        title = NA_character_,
        relation_value = NA_character_ # value controls the vidth of edges so I use relation_value
      )
    } else {
      data.frame(from = entry1, to = entry2, type = type, relation_subtype = xml_attr(
        subnodes,
        "name"
      ), relation_value = xml_attr(subnodes, "value"))
    }
  })

  # Combine into single data frame
  template <- init_empty_edges_df()
  edges_df <- collapse_to_dataframe(template, rels_list)

  message("Parsed ", nrow(edges_df), " relationship edges from KGML file.")

  return(edges_df)
}

#' Download KEGG KGML file for a given pathway ID.
#'
#' @param kgml_file Path to the KGML XML file.  
#' @return A data.frame with the following columns:
#' @export
parse_kgml_reactions <- function(kgml_file) {
  doc <- xml2::read_xml(kgml_file)
  reactions <- xml2::xml_find_all(doc, ".//reaction")

  edges_list <- lapply(reactions, function(reaction) {
    # Reaction attributes
    reaction_id <- xml2::xml_attr(reaction, "id")
    reaction_name <- xml2::xml_attr(reaction, "name")
    reaction_type <- xml2::xml_attr(reaction, "type")

    # Substrates and products
    substrates_nodes <- xml2::xml_find_all(reaction, ".//substrate")
    products_nodes <- xml2::xml_find_all(reaction, ".//product")

    substrates <- if (length(substrates_nodes) > 0) {
      data.frame(
        substrate_id = xml2::xml_attr(substrates_nodes, "id"),
        substrate_name = xml2::xml_attr(substrates_nodes, "name"),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(substrate_id = NA, substrate_name = NA)
    }

    products <- if (length(products_nodes) > 0) {
      data.frame(
        product_id = xml2::xml_attr(products_nodes, "id"),
        product_name = xml2::xml_attr(products_nodes, "name"),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(product_id = NA, product_name = NA)
    }

    # Generate edges: all substrates -> all products
    if (nrow(substrates) > 0 & nrow(products) > 0) {
      edges <- expand.grid( # https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/expand.grid
        from = substrates$substrate_id,
        from_name = substrates$substrate_name,
        to = products$product_id,
        to_name = products$product_name,
        stringsAsFactors = FALSE
      )
      # Add reaction info to edges
      edges <- edges |>
        transform(
          reaction_id = reaction_id,
          reaction_name = reaction_name,
          reaction_type = reaction_type
        )
      return(edges)
    } else {
      return(NULL)
    }
  })

  # Combine all edges into a single data frame
  template <- init_empty_edges_df()
  edges_df <- collapse_to_dataframe(template, edges_list)

  message("Parsed ", nrow(edges_df), " reaction edges from KGML file.")
  return(edges_df)
}

#' Parse KEGG KGML files to extract combined edges data frame.
#'
#' @param file Path to the KGML XML file.
#' @return A data.frame combining both relations and reactions edges.
#' @details
#' This function combines the outputs of `parse_kgml_relations` and
#' `parse_kgml_reactions` to provide a comprehensive edges data frame
#' representing all interactions defined in the KGML file.
#' @export
parse_kgml_edges <- function(file) {
  relations_df <- parse_kgml_relations(file)
  reactions_df <- parse_kgml_reactions(file)

  # Combine both data frames
  edges_df <- rbind(relations_df, reactions_df)
  edges_df <- init_empty_cols_edges(edges_df)

  message("Total edges parsed from KGML file: ", nrow(edges_df))
  return(edges_df)
}

#' Parse KEGG KGML files to extract nodes data frame.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{id}{Unique identifier of the entry within the KGML pathway
#' (from the `id` attribute).}
#'   \item{kegg_name}{The KEGG-specific name or identifier of the entity
#' (from the `name` attribute).
#' This may include one or more KEGG identifiers such as gene IDs,
#' compound IDs, or enzyme EC numbers. Preceded by organism ID}
#'   \item{type}{Type of the node (from the `type` attribute), indicating the biological
#' entity class such as `'gene'`, `'enzyme'`, `'compound'`, `'map'`, `'ortholog'`, or `'group'`.}
#'   \item{link}{URL linking to the KEGG resource for this entry, if available
#' (from the `link` attribute).}
#'   \item{reaction}{Associated reaction ID(s), if any (from the `reaction` attribute).
#' Typically present for enzyme entries.}
#'   \item{graphics_name}{Display name for the entry, taken from the `name`
#'  attribute of the `<graphics>` node.}
#'   \item{label}{Text label for visualization purposes.}
#'   \item{fgcolor}{Foreground color of the graphical element
#'  (from the `fgcolor` attribute of `<graphics>`).}
#'   \item{bgcolor}{Background color of the graphical element
#' (from the `bgcolor` attribute of `<graphics>`).}
#'   \item{graphics_type}{Shape or representation type of the graphical element
#' (from the `type` attribute of `<graphics>`), such as `'rectangle'`, `'circle'`, or `'line'`.}
#'   \item{x}{X-coordinate of the node’s position in the pathway diagram
#' (from the `x` attribute of `<graphics>`).}
#'   \item{y}{Y-coordinate of the node’s position in the pathway diagram
#' (from the `y` attribute of `<graphics>`).}
#'   \item{width}{Width of the graphical element (from the `width` attribute of `<graphics>`).}
#'   \item{height}{Height of the graphical element (from the `height` attribute of `<graphics>`).}
#' }
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr
#'
#' @details
#' The function parses a KEGG KGML (KEGG Markup Language) XML file and extracts all
#' `<entry>` elements, each representing a biological entity in a KEGG pathway diagram.
#' Each node may contain nested `<graphics>` elements defining visual properties
#' (such as position, size, and colors) and `<component>` elements that define group
#' membership for composite entities.
#'
#' The resulting data provides a tidy, one-row-per-entry representation suitable
#' for integration with relational data models or network visualization frameworks
#' (e.g., `igraph`or `visNetwork`).
#'
#' @export
parse_kgml_entries <- function(file) {
  # Read the KGML file
  doc <- xml2::read_xml(file)

  # Find all entry nodes
  entries <- xml2::xml_find_all(doc, ".//entry")

  # Map over each entry (can then have do.call but do.call returns dataframe)
  nodes_list <- lapply(entries, function(entry) {
    graphics_nodes <- xml2::xml_find_all(entry, ".//graphics")
    group_components <- xml2::xml_find_all(entry, ".//component")

    # Base row with initialized attributes
    node_row <- data.frame(
      name = xml2::xml_attr(entry, "id"),
      id = xml2::xml_attr(entry, "id"),
      kegg_name = xml2::xml_attr(entry, "name"),
      type = xml2::xml_attr(entry, "type"),
      link = xml2::xml_attr(entry, "link"),
      reaction = xml2::xml_attr(entry, "reaction"),
      stringsAsFactors = FALSE
    )

    # Extract graphics attributes (only first graphics node is used)
    if (length(graphics_nodes) > 0) {
      if (length(graphics_nodes) > 1) {
        warning(paste("Entry", node_row$id, "has multiple graphics nodes; using the first one."))
      }
      g <- graphics_nodes[[1]]
      node_row$graphics_name <- xml2::xml_attr(g, "name")
      node_row$label <- xml2::xml_attr(g, "name") # for visNetwork
      node_row$fgcolor <- xml2::xml_attr(g, "fgcolor")
      node_row$bgcolor <- xml2::xml_attr(g, "bgcolor")
      node_row$graphics_type <- xml2::xml_attr(g, "type")
      node_row$x <- xml2::xml_attr(g, "x")
      node_row$y <- xml2::xml_attr(g, "y")
      node_row$width <- xml2::xml_attr(g, "width")
      node_row$height <- xml2::xml_attr(g, "height")
    }

    # Extract group components
    if (length(group_components) > 0) {
      node_row$components <- paste(xml2::xml_attr(group_components, "id"),
        collapse = ";"
      )
    }

    node_row
  })

  # Combine all nodes into a single data frame
  template <- init_empty_nodes_df()
  nodes_df <- collapse_to_dataframe(template, nodes_list)
  # Extract KEGG IDs
  nodes_df$KEGG <- vapply(nodes_df$kegg_name, remove_kegg_prefix_str, character(1))
  nodes_df <- init_empty_cols_nodes(nodes_df)
  message("Parsed ", nrow(nodes_df), " nodes from KGML file.")
  return(nodes_df)
}

#' Initialize an empty edges data frame with the required columns.
#' @return An empty data.frame with predefined columns for edges.
#' @noRd
init_empty_edges_df <- function() {
  data.frame(
    from             = character(0),
    to               = character(0),
    type             = character(0),
    relation_subtype = character(0),
    relation_value   = character(0),
    title            = character(0),
    width            = numeric(0),
    color            = character(0),
    arrows           = character(0),
    dashes           = logical(0),
    label            = character(0),
    reaction_id      = character(0),
    reaction_name    = character(0),
    reaction_type    = character(0),
    from_name        = character(0),
    to_name          = character(0),
    stringsAsFactors = FALSE
  )
}

#' Initialize an empty nodes data frame with the required columns.
#' @return An empty data.frame with predefined columns for nodes.
#' @noRd
init_empty_nodes_df <- function() {
  data.frame(
    name             = character(0),
    id               = character(0),
    kegg_name        = character(0),
    type             = character(0),
    link             = character(0),
    reaction         = character(0),
    graphics_name    = character(0),
    label            = character(0),
    fgcolor          = character(0),
    bgcolor          = character(0),
    graphics_type    = character(0),
    x                = character(0),
    y                = character(0),
    choords          = character(0),
    width            = character(0),
    height           = character(0),
    components       = character(0),
    plot_value       = numeric(0),
    source           = character(0),
    color            = character(0),
    text             = character(0),
    group            = character(0),
    fixed            = logical(0),
    widthConstraint  = numeric(0),
    heightConstraint = numeric(0),
    size             = numeric(0),
    shape            = character(0),
    stringsAsFactors = FALSE
  )
}

#' Initialize empty columns for nodes data frame.
#' @param nodes_df Data frame of nodes.
#' @return Data frame with additional empty columns initialized.
#' @noRd
init_empty_cols_nodes <- function(nodes_df) {
  nodes_df$plot_value <- NA_real_
  nodes_df$source <- NA_character_
  nodes_df$color <- nodes_df$bgcolor
  nodes_df$text <- ""
  nodes_df$group <- NA_character_
  nodes_df$fixed <- FALSE
  nodes_df$widthConstraint <- NA_real_
  nodes_df$heightConstraint <- NA_real_
  nodes_df$size <- NA_real_
  nodes_df$shape <- NA_character_

  return(nodes_df)
}

#' Initialize empty columns for edges data frame.
#' @param edges_df Data frame of edges.  
#' @return Data frame with additional empty columns initialized.
#' @noRd
init_empty_cols_edges <- function(edges_df) {
  if (nrow(edges_df) == 0) {
    edges_df$width <- numeric(0)
    edges_df$color <- character(0)
    edges_df$arrows <- character(0)
    edges_df$dashes <- logical(0)
    edges_df$label <- character(0)
  } else {
    edges_df$width <- 1
    edges_df$color <- "gray"
    edges_df$arrows <- "to"
    edges_df$dashes <- FALSE
    edges_df$label <- ""
    edges_df$title <- NA_character_
  }

  return(edges_df)
}

#' Collapse a list of data frames into a single data frame, aligning columns to a template.
#'
#' @param template_df A data frame defining the desired column structure.
#' @param dfs_list A list of data frames to be combined.
#' @return A single data frame combining all input data frames, aligned to the template.
#' @noRd
collapse_to_dataframe <- function(template_df, dfs_list) {
  # Remove NULLs
  dfs_list <- Filter(Negate(is.null), dfs_list)

  # Return empty template if no data frames to combine (or only NULLs)
  if (length(dfs_list) == 0) {
    warning("No entries found in kgml file.")
    return(template_df)
  }

  dfs_aligned <- lapply(dfs_list, function(df) {
    # Add missing columns from template and fill with NA
    missing_cols <- setdiff(names(template_df), names(df))

    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }

    # Reorder columns to match template
    df <- df[, names(template_df), drop = FALSE]
    df
  })
  # Combine all
  do.call(rbind, dfs_aligned)
}
