#' Add tooltips to nodes for visNetwork visualization.
#' @param vertices_df Data frame of nodes with columns: KEGG, label, source, value.
#' @return vertices_df with added 'title' column for tooltips.
#' @details The tooltip includes a button to the specific KEGG entry page.
#' If multiple KEGG IDs are present, they are concatenated with '+' in the URL.
#' It also adds information about the node name, source of differential
#' expression data, and value.
#' @noRd
add_node_tooltip <- function(vertices_df) {
  vertices_df$title <- vapply(
    seq_len(nrow(vertices_df)),
    function(i) {
      switch(tolower(vertices_df$type[i]),
        group = group_node_html(vertices_df[i, , drop = FALSE]),
        ortholog = regular_node_html(vertices_df[i, , drop = FALSE]),
        gene = regular_node_html(vertices_df[i, , drop = FALSE]),
        enzyme = regular_node_html(vertices_df[i, , drop = FALSE]),
        compound = regular_node_html(vertices_df[i, , drop = FALSE]),
        other_node_html(vertices_df[i, , drop = FALSE])
      )
    },
    character(1)
  )

  vertices_df
}


#' Add tooltips to other nodes
#' @param vertices_df Data frame of nodes
#' @return HTML string for the tooltip
#' @noRd
other_node_html <- function(vertices_df) {
  paste0(
    "<h4 style='text-align: center;'>", vertices_df$label, "</h4>",
    "<table>",
    "<tr><th align='left'>Name </th><td>",
    ifelse(is.na(vertices_df$KEGG), "N/A", vertices_df$KEGG),
    "</td></tr>",
    "<tr><th align='left'>ID </th><td>",
    ifelse(is.na(vertices_df$name), "N/A", vertices_df$name),
    "</td></tr>",
    "</table>"
  )
}

#' Add tooltips to group nodes
#' @param vertices_df Data frame of nodes
#' @return HTML string for the tooltip
#' @noRd
group_node_html <- function(vertices_df) {
  paste0(
    "<table>",
    "<tr><th align='left'>Group</th><td>",
    ifelse(
      is.na(vertices_df$group) | vertices_df$group == "",
      "Not part of any group",
      vertices_df$group
    ),
    "</td></tr>",
    "</table>"
  )
}

#' Add tooltips to regular nodes
#' @param vertices_df Data frame of nodes
#' @return HTML string for the tooltip
#' @noRd
regular_node_html <- function(vertices_df) {
  button_html_name <- ifelse(
    is.na(vertices_df$link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", vertices_df$link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#337ab7; border-color:#2e6da4;'>",
      paste0(vertices_df$label, "on KEGG"),
      "</button></a></div>"
    )
  )
  button_html_reaction <- ifelse(
    is.na(vertices_df$reaction_link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", vertices_df$reaction_link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#33b742; border-color:#2e6da4;'>",
      paste0(vertices_df$reaction_name, "on KEGG"),
      "</button></a></div>"
    )
  )
  button_html <- paste0(
    "<div style='display:flex; justify-content:center; gap:6px; margin-top:5px;'>",
    button_html_name,
    button_html_reaction,
    "</div>"
  )
  vertices_df$title <- paste0(
    "<h4 style='text-align: center;'>", vertices_df$label, "</h4>",
    "<table>",
    "<tr><th align='left'>KEGG ID</th><td>",
    ifelse(
      is.na(vertices_df$KEGG),
      "N/A",
      ifelse(
        nchar(vertices_df$KEGG) > 50,
        substr(vertices_df$KEGG, 1, 50),
        vertices_df$KEGG
      )
    ),
    "</td></tr>",
    "<tr><th align='left'>Name</th><td>",
    ifelse(is.na(vertices_df$graphics_name), "N/A", vertices_df$graphics_name),
    "</td></tr>",

    # Source row
    ifelse(
      is.na(vertices_df$de_source) | vertices_df$de_source == "",
      "",
      paste0(
        "<tr><th align='left'>Source</th><td>",
        vertices_df$de_source,
        "</td></tr>"
      )
    ),

    # Value row
    ifelse(
      is.na(vertices_df$de_value) | vertices_df$de_value == "",
      "",
      paste0(
        "<tr><th align='left'>", vertices_df$de_name, "</th><td>",
        round(as.numeric(vertices_df$de_value), 3),
        "</td></tr>"
      )
    ),

    # Group row
    ifelse(
      is.na(vertices_df$group) | vertices_df$group == "",
      "",
      paste0(
        "<tr><th align='left'>Group</th><td>",
        vertices_df$group,
        "</td></tr>"
      )
    ),
    "</table>",
    button_html
  )
}


#' Add tooltips to edges for visNetwork visualization.
#' @param edges_df Data frame of edges with columns: relation_subtype, type, label.
#' @return edges_df with added 'title' column for tooltips.
#' @noRd
add_edge_tooltip <- function(edges_df) {
  base_url <- "https://www.genome.jp/dbget-bin/www_bget?"

  button_html_reaction <- ifelse(
    is.na(edges_df$reaction_name),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", base_url, edges_df$reaction_name, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#337ab7; border-color:#2e6da4;'>",
      "KEGG entry",
      "</button></a></div>"
    )
  )
  # Initialize the title column
  edges_df$title <- ""

  # Relation edges
  idx_relation <- edges_df$type == "relation"
  edges_df$title[idx_relation] <- paste0(
    "<h4 style='text-align: center;'>", edges_df$type[idx_relation], "</h4>",
    "<table>",
    "<tr><th align='left'>Type: </th><td>", edges_df$relation_type[idx_relation], "</td></tr>",
    "<tr><th align='left'>Subtype: </th><td>", edges_df$relation_subtype_name[idx_relation], "</td></tr>",
    "<tr><th align='left'>Label: </th><td>", edges_df$relation_subtype_value[idx_relation], "</td></tr>",
    "</table>"
  )

  # Reaction edges
  idx_reaction <- edges_df$type == "reaction"
  edges_df$title[idx_reaction] <- paste0(
    "<h4 style='text-align: center;'>", edges_df$type[idx_reaction], "</h4>",
    "<table>",
    "<tr><th align='left'>ID: </th><td>", edges_df$reaction_id[idx_reaction], "</td></tr>",
    "<tr><th align='left'>Type: </th><td>", edges_df$reaction_type[idx_reaction], "</td></tr>",
    "<tr><th align='left'>Name: </th><td>", edges_df$reaction_name[idx_reaction], "</td></tr>",
    "</table>",
    button_html_reaction[idx_reaction] 
  )

  # Line edges
  idx_line <- edges_df$type == "line"
  edges_df$title[idx_line] <- paste0(
    "<h4 style='text-align: center;'>", edges_df$type[idx_line], "</h4>",
    "<table>",
    "<tr><th align='left'>Name: </th><td>", edges_df$reaction_name[idx_line], "</td></tr>",
    "</table>",
    button_html_reaction[idx_line] 
  )

  return(edges_df)
}

