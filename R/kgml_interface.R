#' Parse KEGG KGML files to extract relations and edges data frames.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{from}{The ID of the source entry (node) in the pathway, corresponding to the `entry1` attribute in the KGML `<relation>` tag.}
#'   \item{to}{The ID of the target entry (node) in the pathway, corresponding to the `entry2` attribute in the KGML `<relation>` tag.}
#'   \item{type}{The general type of relationship between the two entries (e.g., `"ECrel"`, `"PPrel"`, `"GErel"`, `"PCrel"`, `"maplink"`). These types describe the biological nature of the connection, such as enzyme-enzyme relation or protein-protein interaction.}
#'   \item{subtype}{A more specific subtype of the relation, derived from the `<subtype>` child elements of the `<relation>` node (e.g., `"activation"`, `"inhibition"`, `"expression"`, `"compound"`). If no subtype is defined, this will be `NA`.}
#'   \item{rel_value}{A categorical value associated with the subtype, which encodes the apearence of the arrow (from the `value` attribute of `<subtype>`). If not present, this will be `NA`.}
#' }
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr
#'
#' @details
#' The function reads a KEGG KGML (KEGG Markup Language) file, which encodes pathway
#' information as XML, and extracts all `<relation>` elements that describe the
#' interactions or relationships between entities in the pathway. Each `<relation>` may
#' contain one or more `<subtype>` elements that provide additional details about
#' the interaction. The result is a tidy data.frame suitable for network analysis or
#' visualization, where each row represents one relation–subtype pair.
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
        subtype = NA_character_,
        rel_value = NA_character_ # value controls the vidth of edges
      )
    } else {
      data.frame(
        from = entry1,
        to = entry2,
        type = type,
        subtype = xml_attr(subnodes, "name"),
        rel_value = xml_attr(subnodes, "value")
      )
    }
  })
  edges_df <- do.call(rbind, rels_list)

  # Default to empty data.frame if no relations found
  if (is.null(edges_df)) {
    warning("No relations found in KGML file.")
    edges_df <- data.frame(
      from = character(0),
      to = character(0),
      type = character(0),
      subtype = character(0),
      rel_value = character(0)
    )
  }

  message("Parsed ", nrow(edges_df), " edges from KGML file.")

  return(edges_df)
}


#' Parse KEGG KGML files to extract nodes data frame.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{id}{Unique identifier of the entry within the KGML pathway (from the `id` attribute).}
#'   \item{kegg_name}{The KEGG-specific name or identifier of the entity (from the `name` attribute). This may include one or more KEGG identifiers such as gene IDs, compound IDs, or enzyme EC numbers. Preceded by organism ID}
#'   \item{type}{Type of the node (from the `type` attribute), indicating the biological entity class such as `"gene"`, `"enzyme"`, `"compound"`, `"map"`, `"ortholog"`, or `"group"`.}
#'   \item{link}{URL linking to the KEGG resource for this entry, if available (from the `link` attribute).}
#'   \item{reaction}{Associated reaction ID(s), if any (from the `reaction` attribute). Typically present for enzyme entries.}
#'   \item{graphics_name}{Display name for the entry, taken from the `name` attribute of the `<graphics>` node.}
#'   \item{label}{Text label for visualization purposes.}
#'   \item{fgcolor}{Foreground color of the graphical element (from the `fgcolor` attribute of `<graphics>`).}
#'   \item{bgcolor}{Background color of the graphical element (from the `bgcolor` attribute of `<graphics>`).}
#'   \item{graphics_type}{Shape or representation type of the graphical element (from the `type` attribute of `<graphics>`), such as `"rectangle"`, `"circle"`, or `"line"`.}
#'   \item{x}{X-coordinate of the node’s position in the pathway diagram (from the `x` attribute of `<graphics>`).}
#'   \item{y}{Y-coordinate of the node’s position in the pathway diagram (from the `y` attribute of `<graphics>`).}
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
      name = xml2::xml_attr(entry, "id"), # set name because thats how igraph builds connections
      id = xml2::xml_attr(entry, "id"),
      kegg_name = xml2::xml_attr(entry, "name"),
      type = xml2::xml_attr(entry, "type"),
      link = xml2::xml_attr(entry, "link"),
      reaction = xml2::xml_attr(entry, "reaction"),
      graphics_name = NA_character_,
      label = NA_character_,
      fgcolor = NA_character_,
      bgcolor = NA_character_,
      graphics_type = NA_character_,
      x = NA_character_,
      y = NA_character_,
      width = NA_character_,
      height = NA_character_,
      components = NA_character_
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
      node_row$components <- paste(xml2::xml_attr(group_components, "id"), collapse = ";")
    }

    node_row
  })
  nodes_df <- do.call(rbind, nodes_list)
  message("Parsed ", nrow(nodes_df), " nodes from KGML file.")
  return(nodes_df)
}

#' Add columns to nodes_df with default values
#' @param nodes_df Data frame of nodes parsed from KGML
#' @return Data frame with additional columns for styling and DE results
add_columns_nodes_df <- function(nodes_df) {
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
