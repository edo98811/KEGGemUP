# https://rpubs.com/huanfaChen/261838

#' @noRd
kegg_node_defaults <- function() {
  list(
    # Core KGML attributes
    name = NA_character_,
    kegg_entry_id = NA_character_,
    KEGG_no_prefix = "",
    type = NA_character_,
    link = NA_character_,
    reaction = NA_character_,
    reaction_link = NA_character_,
    reaction_label = NA_character_,
    graphics_name = NA_character_,
    node_name = NA_character_,
    KEGG = NA_character_,
    components = NA_character_,
    line_id = NA_character_,
    coords = NA_character_,
    graphics_type = NA_character_,
    ids_for_mapping = "",

    # Visualization attributes
    label = "", # minimum necessary
    fgcolor = NA_character_,
    bgcolor = NA_character_,
    x = NA_real_,
    y = NA_real_,
    width = NA_real_,
    height = NA_real_,
    source = NA_character_,
    fgcolor = NA_character_,
    bgcolor = NA_character_,
    group = NA_character_,
  )
}

#' @noRd
kegg_edge_defaults <- function() {
  list(
    # Core KGML attributes
    from = NA_character_,
    to = NA_character_,
    type = NA_character_,
    relation_type = NA_character_,
    relation_subtype_name = NA_character_,
    relation_subtype_value = NA_character_,
    reaction_alt_name_substrate = NA_character_,
    reaction_alt_name_product = NA_character_,
    reaction_id = NA_character_,
    reaction_name = NA_character_,
    reaction_type = NA_character_,
    reaction_from_name = NA_character_,
    reaction_to_name = NA_character_,
    point_index = NA_integer_,
  )
}


#' @noRd
node_defaults <- function() {
  list(
    # Basic attributes
    name = NA_integer_, # kegg unique node id
    label = NA_character_,
    title = NA_character_,
    group = NA_character_,
    link = NA_character_,
    de_value = NA_real_,
    de_source = NA_character_,
    feature_id_1 = NA_character_, # KEGG: KEGG ids
    feature_id_2 = NA_character_, # KEGG: graphics_name
    original_shape = NA_character_,
    ids_for_mapping = "", # KEGG: ENTREZID

    # Style attributes
    color = "white",
    size = 25,
    x = NA_integer_,
    y = NA_integer_,
    fixed = FALSE,
    width = NA_integer_, # for visNetwork
    height = NA_integer_, # for visNetwork
    borderRadius = NA_integer_, # for visNetwork
    shape = "vrectangle",
    text = ""
  )
}

#' @noRd
edge_defaults <- function() {
  list(
    # Basic attributes
    from = NA_character_,
    to = NA_character_,
    id = uuid::UUIDgenerate(),
    type = NA_character_,
    link = NA_character_,
    name_1 = NA_character_, # KEGG: Name
    name_2 = NA_character_, # KEGG: Type
    name_3 = NA_character_, # KEGG: Subtype

    # Style attributes
    directed = TRUE,
    color = "gray",
    width = 1,
    value = NA_real_,
    label = "",
    lty = "solid",
    arrows = "to",
    dashes = FALSE,
    title = NA_character_
  )
}

#' @noRd
kegg_to_general_node_map <- function() {
  c(
    # Identity / annotation
    name = "name",
    label = "label",
    link = "link",
    type = "type",
    KEGG = "feature_id_1", # KEGG  ID
    graphics_name = "feature_id_2", # the long name used in KEGG graphics
    KEGG_no_prefix = "ids_for_mapping",
    original_shape = "graphics_type",

    # Layout / geometry
    x = "x",
    y = "y",
    width = "width",
    height = "height",

    # Styling
    text = "text"
  )
}


#' @noRd
kegg_to_general_edge_map <- function() {
  c(
    # Topology
    from = "from",
    to = "to",

    # Identity / annotation
    type = "type", # KEGG: Type (e.g., PPrel, GErel)
    link = "link",
    relation_type = "name_1",
    relation_subtype_name = "name_2",
    relation_subtype_value = "name_3",
    reaction_alt_name_substrate = "",
    reaction_alt_name_product = "",
    reaction_id = "name_1",
    reaction_name = "name_2",
    reaction_type = "name_3",
  )
}
