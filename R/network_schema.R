# https://rpubs.com/huanfaChen/261838

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
    feature_id_1 = NA_character_, # KEGG: KEGG IDs
    feature_id_2 = NA_character_, # KEGG: ENTREZID
    feature_id_3 = NA_character_, # KEGG: graphics_name

    # Style attributes
    color = NA_character_,
    size = NA_integer_,
    x = NA_integer_,
    y = NA_integer_,
    fixed = FALSE,
    widthConstraint = NA_integer_, # for visNetwork
    heightConstraint = NA_integer_, # for visNetwork
    borderRadius = NA_integer_, # for visNetwork
    shape = NA_character_,
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
    subtype = NA_character_,
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
kegg_node_defaults <- function() {
  list(
    # Core KGML attributes
    name = NA_character_,
    id = NA_character_,
    kegg_name = NA_character_,
    type = NA_character_,
    link = NA_character_,
    reaction = NA_character_,
    graphics_name = NA_character_,
    KEGG = NA_character_,
    components = NA_character_,

    # Visualization attributes
    label = NA_character_,
    fgcolor = NA_character_,
    bgcolor = NA_character_,
    graphics_type = NA_character_,
    x = NA_real_,
    y = NA_real_,
    coords = NA_character_,
    width = NA_real_,
    height = NA_real_,
    plot_value = NA_real_,
    source = NA_character_,
    color = NA_character_,
    borderRadius = NA_integer_, # for visNetwork
    text = "",
    group = NA_character_,
    fixed = FALSE,
    widthConstraint = NA_real_,
    heightConstraint = NA_real_,
    size = NA_real_,
    shape = NA_character_
  )
}

#' @noRd
kegg_edge_defaults <- function() {
  list(
    # Core KGML attributes
    from = NA_character_,
    to = NA_character_,
    type = NA_character_,
    relation_subtype = NA_character_,
    relation_value = NA_character_,
    reaction_id = NA_character_,
    reaction_name = NA_character_,
    reaction_type = NA_character_,
    from_name = NA_character_,
    to_name = NA_character_,

    # Visualization attributes
    title = NA_character_,
    width = 1,
    color = "gray",
    arrows = "to",
    dashes = FALSE,
    label = ""
  )
}

#' @noRd
kegg_to_general_node_map <- function() {
  c(
    # Identity / annotation
    name           = "name",
    label          = "label",
    link           = "link",
    KEGG           = "feature_id_2", # typically ENTREZID
    id             = "feature_id_1", # typically KEGG gene ID
    graphics_name  = "feature_id_3",

    # Layout / geometry
    x              = "x",
    y              = "y",
    width          = "widthConstraint",
    height         = "heightConstraint",
    borderRadius   = "borderRadius",
    size           = "size",
    shape          = "shape",
    fixed          = "fixed",

    # Styling
    text           = "text"
  )
}


#' @noRd
kegg_to_general_edge_map <- function() {
  c(
    # Topology
    entry1   = "from",
    entry2   = "to",

    # Identity / annotation
    name     = "name_1", # KEGG: Name
    type     = "name_2", # KEGG: Type (e.g., PPrel, GErel)
    subtype  = "name_3", # KEGG: Subtype (e.g., activation)
    link     = "link",
    id       = "id",

    # Semantics
    type     = "type",
    subtype  = "subtype",

    # Visualization
    directed = "directed",
    color    = "color",
    width    = "width",
    value    = "value",
    label    = "label",
    lty      = "lty",
    arrows   = "arrows",
    dashes   = "dashes",
    title    = "title"
  )
}
