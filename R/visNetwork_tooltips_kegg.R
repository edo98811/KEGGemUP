#' Add tooltips to nodes for visNetwork visualization.
#' @param nodes_df Data frame of nodes with columns: KEGG, label, source, value.
#' @return nodes_df with added 'title' column for tooltips.
#' @details The tooltip includes a button to the specific KEGG entry page.
#' If multiple KEGG IDs are present, they are concatenated with '+' in the URL.
#' It also adds information about the node name, source of differential
#' expression data, and value.
#' @noRd
add_node_tooltip <- function(nodes_df) {
  button_html_name <- ifelse(
    is.na(nodes_df$link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", nodes_df$link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#337ab7; border-color:#2e6da4;'>",
      "KEGG entry",
      "</button></a></div>"
    )
  )
  button_html_reaction <- ifelse(
    is.na(nodes_df$reaction_link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", nodes_df$reaction_link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:##33b742; border-color:#2e6da4;'>",
      "KEGG entry",
      "</button></a></div>"
    )
  )
  button_html <- paste0(
    "<div style='display:flex; justify-content:center; gap:6px; margin-top:5px;'>",
    button_html_name,
    button_html_reaction,
    "</div>"
  )
  nodes_df$title <- ifelse(
    nodes_df$kegg_name == "undefined",

    # Group placeholder node
    paste0(
      "<table>",
      "<tr><th align='left'>Group</th><td>",
      ifelse(
        is.na(nodes_df$group) | nodes_df$group == "",
        "Not part of any group",
        nodes_df$group
      ),
      "</td></tr>",
      "</table>"
    ),

    # Regular node
    paste0(
      "<h4 style='text-align: center;'>", nodes_df$label, "</h3>",
      "<table>",
      "<tr><th align='left'>KEGG ID  </th><td>",
      ifelse(
        nchar(nodes_df$kegg_name) > 50,
        substr(nodes_df$kegg_name, 1, 50),
        nodes_df$kegg_name
      ),
      "</td></tr>",
      "<tr><th align='left'>Name</th><td>",
      ifelse(is.na(nodes_df$feature_id_2), "N/A", nodes_df$feature_id_2),
      "</td></tr>",
      "<tr><th align='left'>Source</th><td>",
      ifelse(is.na(nodes_df$de_source), "N/A", nodes_df$de_source),
      "</td></tr>",
      "<tr><th align='left'>Value</th><td>",
      ifelse(
        is.na(nodes_df$plot_value), "",
        format(round(as.numeric(nodes_df$plot_value), 3), nsmall = 3)
      ),
      "</td></tr>",
      "<tr><th align='left'>Group</th><td>",
      ifelse(
        is.na(nodes_df$group) | nodes_df$group == "",
        "Not belonging to any group",
        nodes_df$group
      ),
      "</td></tr>",
      "</table>",
      button_html
    )
  )


  return(nodes_df)
}

#' Add tooltips to edges for visNetwork visualization.
#' @param edges_df Data frame of edges with columns: relation_subtype, type, label.
#' @return edges_df with added 'title' column for tooltips.
#' @noRd
add_edge_tooltip <- function(edges_df) {
  edges_df$title <- paste0(
    "<table>",
    "<tr><th align='left'>Type</th><td>", edges_df$type, "</td></tr>",
    "<tr><th align='left'>",
    ifelse(edges_df$type == "reaction", "ID", "Relation Type"),
    "</th><td>",
    ifelse(edges_df$type == "reaction", edges_df$ID, edges_df$relation_subtype),
    "</td></tr>",
    "<tr><th align='left'>Name</th><td>",
    ifelse(is.na(edges_df$name) | edges_df$name == "", "N/A", edges_df$name),
    "</td></tr>",
    "<tr><th align='left'>Subtype</th><td>",
    ifelse(is.na(edges_df$subtype) | edges_df$subtype == "", "N/A", edges_df$subtype),
    "</td></tr>",
    "</table>"
  )
  return(edges_df)
}
