#' Parse KEGG KGML files to extract relations and edges data frames.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{from}{The ID of the source entry (node) in the pathway, corresponding to the `entry1` attribute in the KGML `<relation>` tag.}
#'   \item{to}{The ID of the target entry (node) in the pathway, corresponding to the `entry2` attribute in the KGML `<relation>` tag.}
#'   \item{type}{The general type of relationship between the two entries (e.g., `"ECrel"`, `"PPrel"`, `"GErel"`, `"PCrel"`, `"maplink"`). These types describe the biological nature of the connection, such as enzyme-enzyme relation or protein-protein interaction.}
#'   \item{subtype}{A more specific subtype of the relation, derived from the `<subtype>` child elements of the `<relation>` node (e.g., `"activation"`, `"inhibition"`, `"expression"`, `"compound"`). If no subtype is defined, this will be `NA`.}
#'   \item{rel_value}{A categorical value associated with the subtype, which encodes the apearence of the arrow (from the `value` attribute of `<subtype>`). If not present, this will be `NA`.}
#' }
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
#'
#' @details
#' The function reads a KEGG KGML (KEGG Markup Language) file, which encodes pathway
#' information as XML, and extracts all `<relation>` elements that describe the 
#' interactions or relationships between entities in the pathway. Each `<relation>` may
#' contain one or more `<subtype>` elements that provide additional details about
#' the interaction. The result is a tidy tibble suitable for network analysis or
#' visualization, where each row represents one relation–subtype pair.
#' @examples
#' \dontrun{
#' relations_df <- parse_kgml_relations("pathway.xml")
#' head(relations_df)
#' }
#' @export
parse_kgml_relations <- function(file) {
  doc <- read_xml(file)

  rels <- xml_find_all(doc, ".//relation")

  purrr::map_dfr(rels, function(rel) {
    entry1 <- xml_attr(rel, "entry1")
    entry2 <- xml_attr(rel, "entry2")
    type <- xml_attr(rel, "type")
    subnodes <- xml_find_all(rel, ".//subtype")

    if (length(subnodes) == 0) {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        subtype = NA_character_,
        rel_value = NA_character_ # value controls the vidth of edges
      )
    } else {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        subtype = xml_attr(subnodes, "name"),
        rel_value = xml_attr(subnodes, "value")
      )
    }
  })
}


#' Parse KEGG KGML files to extract nodes data frame.
#'
#' @param file Path to the KGML XML file.
#'
#' @return A tibble with the following columns:
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
#'   \item{value}{Will be filled up later with numerical values mapped to the node (initialized as `NA`).}
#'   \item{source}{Will be filled up later with text indicating the source of the numerical values mapped to the node (initialized as `NA`).}
#'   \item{color}{Convenience column duplicating `bgcolor`, used for visualization with `visNetwork` or `igraph`.}
#'   \item{text}{text field (initialized as empty string) for additional infos on mapping.}
#'   \item{components}{Semicolon-separated list of component IDs for entries of type `"group"`, extracted from nested `<component>` nodes.}
#'   \item{group}{Reserved column for grouping information (initialized as NA).}
#' }
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
#'
#' @details
#' The function parses a KEGG KGML (KEGG Markup Language) XML file and extracts all
#' `<entry>` elements, each representing a biological entity in a KEGG pathway diagram.
#' Each node may contain nested `<graphics>` elements defining visual properties 
#' (such as position, size, and colors) and `<component>` elements that define group
#' membership for composite entities.
#'
#' The resulting tibble provides a tidy, one-row-per-entry representation suitable
#' for integration with relational data models or network visualization frameworks
#' (e.g., `igraph`or `visNetwork`).
#'
#' @examples
#' \dontrun{
#' nodes_df <- parse_kgml_entries("pathway.xml")
#' head(nodes_df)
#' }
#' @export
parse_kgml_entries <- function(file) {
  # Read the KGML file
  doc <- xml2::read_xml(file)

  # Find all entry nodes
  entries <- xml2::xml_find_all(doc, ".//entry")

  # Map over each entry (can then have do.call but do.call returns dataframe)
  purrr::map_dfr(entries, function(entry) {
    graphics_nodes <- xml2::xml_find_all(entry, ".//graphics")
    group_components <- xml2::xml_find_all(entry, ".//component")

    # Base row with initialized attributes
    node_row <- tibble::tibble(
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
      value = NA_real_,
      source = NA_character_,
      color = NA_character_,
      text = "",
      components = NA_character_,
      group = NA_character_
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
      node_row$color <- xml2::xml_attr(g, "bgcolor")
    }

    # Extract group components
    if (length(group_components) > 0) {
      node_row$components <- paste(xml2::xml_attr(group_components, "id"), collapse = ";")
    }

    node_row
  })
}

# #' Parse KEGG KGML files to extract nodes and edges data frames.
# #'
# #' @param file Path to the KGML XML file.
# #' @return A tibble with columns: id, name, type, link, reaction, graphics_name, fgcolor, bgcolor, graphics_type, x, y, width, height.
# #' @importFrom xml2 read_xml xml_find_all xml_attr
# #' @importFrom tibble tibble
# #' @importFrom purrr map_dfr
# parse_kgml_entries <- function(file) {
#   # read the KGML file
#   doc <- read_xml(file)

#   # find all entry nodes
#   entries <- xml_find_all(doc, ".//entry")

#   # map over each entry
#   purrr::map_dfr(entries, function(entry) {
#     entry_id <- xml_attr(entry, "id")
#     entry_name <- xml_attr(entry, "name")
#     entry_type <- xml_attr(entry, "type")
#     entry_link <- xml_attr(entry, "link")
#     entry_reaction <- xml_attr(entry, "reaction")

#     graphics_nodes <- xml_find_all(entry, ".//graphics")
#     group_components <- xml_find_all(entry, ".//component")

#     # if no graphics, create a single row with NAs
#     if (length(graphics_nodes) == 0) {
#       tibble(
#         name = entry_id,
#         kegg_name = entry_name,
#         type = entry_type,
#         link = entry_link,
#         reaction = entry_reaction,
#         graphics_name = NA_character_,
#         label = NA_character_,
#         fgcolor = NA_character_,
#         bgcolor = NA_character_,
#         graphics_type = NA_character_,
#         x = NA_character_,
#         y = NA_character_,
#         width = NA_character_,
#         height = NA_character_,
#         value = NA_real_,
#         source = NA_character_,
#         color = NA_character_,
#         text = "",
#         components = NA_character_
#       )
#     } else {
#       purrr::map_dfr(graphics_nodes, function(g) {
#         tibble(
#           name = entry_id,
#           kegg_name = entry_name,
#           type = entry_type,
#           link = entry_link,
#           reaction = entry_reaction,
#           graphics_name = xml_attr(g, "name"),
#           label = xml_attr(g, "name"), # for visNetwork
#           fgcolor = xml_attr(g, "fgcolor"),
#           bgcolor = xml_attr(g, "bgcolor"),
#           graphics_type = xml_attr(g, "type"),
#           x = xml_attr(g, "x"),
#           y = xml_attr(g, "y"),
#           width = xml_attr(g, "width"),
#           height = xml_attr(g, "height"),
#           value = NA_real_,
#           source = NA_character_,
#           color = xml_attr(g, "bgcolor"),
#           text = ""
#           components = paste(xml_attr(entry_component, "id"), collapse = ";")
#         )
#       })
#     }
#   })
# }


# parse_kgml_entries <- function(file) {
#   # read the KGML file
#   doc <- read_xml(file)

#   # find all entry nodes
#   entries <- xml_find_all(doc, ".//entry")
#   graphics <- xml_find_all(entries, ".//graphics")

#   # extract attributes into a data.frame
#   df <- tibble(
#     id = xml_attr(entries, "id"),
#     name = xml_attr(entries, "name"),
#     type = xml_attr(entries, "type"),
#     link = xml_attr(entries, "link"),
#     reaction = xml_attr(entries, "reaction"),
#     graphics_name = xml_attr(graphics, "name"),
#     label = xml_attr(graphics, "name"), # for visNetwork
#     fgcolor = xml_attr(graphics, "fgcolor"),
#     bgcolor = xml_attr(graphics, "bgcolor"),
#     graphics_type = xml_attr(graphics, "type"),
#     x = xml_attr(graphics, "x"),
#     y = xml_attr(graphics, "y"),
#     width = xml_attr(graphics, "width"),
#     height = xml_attr(graphics, "height"),
#     value = NA,
#     source = NA_character_,
#     color = xml_attr(graphics, "bgcolor"),
#     text = NA_character_
#   )
#   # Example entry:
#   #     <entry id="184" name="ko:K15359 ko:K18276" type="ortholog" reaction="rn:R09472"
#   # link="https://www.kegg.jp/dbget-bin/www_bget?K15359+K18276">
#   # <graphics name="K15359..." fgcolor="#000000" bgcolor="#FFFFFF"
#   #      type="rectangle" x="303" y="561" width="46" height="17"/>

#   return(df)
# }
