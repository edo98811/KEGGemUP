# https://rpubs.com/huanfaChen/261838

#' Default KEGG node attributes
#' @return A list of default KEGG node attributes
#' @noRd
kegg_node_defaults <- function() {
  list(
    # Core KGML attributes
    name = NA_character_,
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
    point_index = NA_integer_,

    # Data mapping attributes
    de_value = NA_real_,
    de_text = NA_character_,
    de_source = NA_character_,
    ids_for_mapping = "",

    # Visualization attributes
    label = "",
    x = NA_integer_,
    y = NA_integer_,
    width = 1,
    height = 1,
    source = NA_character_,
    fgcolor = NA_character_,
    bgcolor = NA_character_,
    group = NA_character_,
    vertex.color = "white",
    size = 25,
    fixed = TRUE,
    shape = "vrectangle",
    text = ""
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
    # Style attributes
    directed = TRUE,
    color = "gray",
    width = 1,
    value = NA_real_,
    label = "",
    lty = 1,
    arrow.mode = "to",
    dashes = FALSE,
    title = NA_character_
  )
}

