#' Add group information to nodes based on 'undefined' groups.
#' @param nodes_df Data frame of nodes with columns: id, kegg_name, components.
#' @return nodes_df with updated 'group' column.
#' @noRd
add_group <- function(nodes_df) {
  # The nodes that have as kegg name 'undefined' are group nodes
  undefined_idx <- which(!is.na(nodes_df$kegg_name) & nodes_df$kegg_name == "undefined")

  # If there are no undefined nodes, return original df
  if (length(undefined_idx) == 0) {
    return(nodes_df)
  }

  # Subset undefined nodes
  undefined_nodes <- nodes_df[undefined_idx, , drop = FALSE]

  # Iterate over undefined nodes using indeces
  for (i in seq_len(nrow(undefined_nodes))) {
    # If no components in group (empty), skip
    if (is.na(undefined_nodes$components[i]) || undefined_nodes$components[i] ==
      "") {
      next
    } # If group is NA

    # Get component ids and add the group label to them, add the group node
    # itself to this list
    ids <- strsplit(undefined_nodes$components[i], ";", fixed = TRUE)[[1]]
    ids <- append(ids, undefined_nodes$id[i])

    # Make group label group_label <- paste0('group_',
    # undefined_nodes$id[i])

    group_elements <- nodes_df$label[nodes_df$id %in% ids]
    group_label <- paste(group_elements[1:length(group_elements) - 1], collapse = ";")

    # Assign group label to nodes_df
    nodes_df[nodes_df$id %in% ids, "group"] <- group_label
  }

  return(nodes_df)
}


#' Add tooltips to nodes for visNetwork visualization.
#' @param nodes_df Data frame of nodes with columns: KEGG, label, source, value.
#' @return nodes_df with added 'title' column for tooltips.
#' @details The tooltip includes a button to the specific KEGG entry page.
#' If multiple KEGG IDs are present, they are concatenated with '+' in the URL.
#' It also adds information about the node name, source of differential
#' expression data, and value.
#' @noRd
add_tooltip <- function(nodes_df) {
  button_html <- ifelse(
    is.na(nodes_df$kegg_name) | is.na(nodes_df$link) | nodes_df$kegg_name == "",
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", nodes_df$link, "' target='_blank'>",
      "<button type='button' style='color:#fff; background-color:#337ab7; border-color:#2e6da4;'>",
      "KEGG entry",
      "</button></a></div>"
    )
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
      ifelse(is.na(nodes_df$graphics_name), "N/A", nodes_df$graphics_name),
      "</td></tr>",
      "<tr><th align='left'>Source</th><td>",
      ifelse(is.na(nodes_df$source), "N/A", nodes_df$source),
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
    "Type: ", edges_df$relation_type, "<br>",
    "relation_subtype: ", edges_df$relation_subtype, "<br>",
    "Label: ", ifelse(edges_df$label == "", "N/A", edges_df$label)
  )
  return(edges_df)
}

#' Style nodes based on their type for visNetwork visualization.
#' @param nodes_df Data frame of nodes with a column 'type'.
#' @param node_size_multiplier Numeric factor to scale node sizes (default: 1.2).
#' @return nodes_df with added visual styling columns:
#' shape, fixed, widthConstraint, heightConstraint, size.
#' @noRd
style_nodes <- function(nodes_df, node_size_multiplier = 1.2) {
  # Base visual settings (for vinetwork)
  nodes_df$shape[nodes_df$type == "compound"] <- "dot"
  nodes_df$shape[nodes_df$type != "compound"] <- "box"
  # Set size constraints for non-compound nodes (compute numeric vectors
  # first)
  widths_num <- as.numeric(nodes_df$width) * node_size_multiplier
  heights_num <- as.numeric(nodes_df$height) * node_size_multiplier
  non_comp_idx <- which(!is.na(nodes_df$type) & nodes_df$type != "compound")

  if (length(non_comp_idx) > 0) {
    nodes_df$widthConstraint[non_comp_idx] <- widths_num[non_comp_idx]
    nodes_df$heightConstraint[non_comp_idx] <- heights_num[non_comp_idx]
  }

  # Group nodes set dimension to one (very small)
  undef_idx <- which(!is.na(nodes_df$kegg_name) & nodes_df$kegg_name == "undefined")
  if (length(undef_idx) > 0) {
    nodes_df$widthConstraint[undef_idx] <- 1
    nodes_df$heightConstraint[undef_idx] <- 1
  }
  # Dot nodes size
  dot_idx <- which(nodes_df$shape == "dot")
  if (length(dot_idx) > 0) {
    nodes_df$size[dot_idx] <- 7
  }
  # Line point nodes
  line_point_idx <- which(nodes_df$type == "line_point")
  if (length(line_point_idx) > 0) {
    nodes_df$shape[line_point_idx] <- "dot"
    nodes_df$size[line_point_idx] <- 1
    nodes_df$color[line_point_idx] <- "transparent"
    nodes_df$fixed[line_point_idx] <- TRUE
  }

  return(nodes_df)
}


#' Style edges based on their relation_subtype for visNetwork visualization.
#' @param edges_df Data frame of edges with a column 'relation_subtype'.
#' @return edges_df with added visual styling columns: color, dashes, arrows, label.
#' @noRd
style_edges <- function(edges_df) {
  # possible relation_subtypes and their styles
  # name	value	ECrel	PPrel	GErel	Explanation
  # compound	Entry element id attribute value for compound.	*	*		shared with two successive reactions (ECrel) or intermediate of two interacting proteins (PPrel)
  # hidden compound	Entry element id attribute value for hidden compound.	*			shared with two successive reactions but not displayed in the pathway map
  # activation	-->		*		positive and negative effects which may be associated with molecular information below
  # inhibition	--|		*
  # expression	-->			*	interactions via DNA binding
  # repression	--|			*
  # indirect effect	..>		*	*	indirect effect without molecular details
  # state change	...		*		state transition
  # binding/association	---		*		association and dissociation
  # dissociation	-+-		*
  # missing interaction	-/-		*	*	missing interaction due to mutation, etc.
  # phosphorylation	+p		*		molecular events
  # dephosphorylation	-p		*
  # glycosylation	+g		*
  # ubiquitination	+u		*
  # methylation	+m		*

  edge_style_map <- list(
    compound = list(color = "black", dashes = FALSE, arrows = "to", label = ""),
    hidden_compound = list(color = "lightgray", dashes = FALSE, arrows = "to", label = ""),
    activation = list(color = "red", dashes = FALSE, arrows = "to", label = ""),
    inhibition = list(color = "blue", dashes = FALSE, arrows = "tee", label = ""),
    expression = list(color = "red", dashes = TRUE, arrows = "to", label = ""),
    repression = list(color = "blue", dashes = TRUE, arrows = "tee", label = ""),
    indirect_effect = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    state_change = list(color = "gray", dashes = TRUE, arrows = "", label = ""),
    binding_association = list(color = "black", dashes = TRUE, arrows = "", label = ""),
    dissociation = list(color = "gray", dashes = TRUE, arrows = "to", label = ""),
    missing_interaction = list(color = "gray", dashes = TRUE, arrows = "to", label = "-/-"),
    phosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+p"),
    dephosphorylation = list(color = "black", dashes = FALSE, arrows = "to", label = "-p"),
    glycosylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+g"),
    ubiquitination = list(color = "black", dashes = FALSE, arrows = "to", label = "+u"),
    methylation = list(color = "black", dashes = FALSE, arrows = "to", label = "+m"),
    others_unknown = list(color = "black", dashes = TRUE, arrows = "to", label = "?"),
    group_relation = list(color = "transparent", dashes = TRUE, arrows = "", label = ""),
    # For reactions
    reversible = list(color = "black", dashes = TRUE, arrows = "", label = ""),
    irreversible = list(color = "black", dashes = TRUE, arrows = "", label = ""),
    line = list(color = "black", dashes = FALSE, arrows = "", label = "")
  )

  # https://builtin.com/data-science/and-in-r#:~:text=The%20single%20sign%20version%20%7C%20returns,first%20element%20of%20each%20vector.
  edges_df$relation_subtype <- tolower(edges_df$relation_subtype)
  edges_df$relation_subtype <- gsub("[/ ]", "_", edges_df$relation_subtype)
  edges_df$relation_subtype[is.na(edges_df$relation_subtype) |
    !(edges_df$relation_subtype %in% names(edge_style_map))] <- "others_unknown"

  # Vectorized assignment
  edges_df$color <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$color, character(1))
  edges_df$dashes <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$dashes, logical(1))
  edges_df$arrows <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$arrows, character(1))
  edges_df$label <- vapply(edges_df$relation_subtype, function(x) edge_style_map[[x]]$label, character(1))


  return(edges_df)
}


#' Scale node dimensions for better visualization.
#' @param nodes_df Data frame of nodes with x and y coordinates.
#' @param factor Scaling factor (default: 2).
#'
#' @return nodes_df with scaled x and y coordinates.
#' @noRd
scale_dimensions <- function(nodes_df, factor = 2) {
  # Scale x and y coordinates to make the graph look nicer
  nodes_df$x <- as.numeric(nodes_df$x) * factor
  nodes_df$y <- as.numeric(nodes_df$y) * factor

  return(nodes_df)
}

#' Add compound names to compound nodes in the nodes data frame.
#' @param nodes_df Data frame of nodes with a column 'type' indicating node type.
#' @param bfc BiocFileCache object for caching KEGG compound mappings.
#' @return Updated nodes data frame with compound names added to compound nodes.
#' @importFrom BiocFileCache BiocFileCache
#' @noRd
add_compound_names <- function(nodes_df, bfc) {
  idx <- which(!is.na(nodes_df$type) & nodes_df$type == "compound")

  if (length(idx) == 0) {
    return(nodes_df)
  }

  compounds_in_graph <- as.character(nodes_df$KEGG)
  compounds_in_graph[is.na(compounds_in_graph)] <- ""
  compounds_in_graph <- compounds_in_graph[idx]

  compounds <- get_kegg_db(db_name = "compound", bfc = bfc) # expect named vector mapping KEGG id -> name
  glycan <- get_kegg_db(db_name = "glycan", bfc = bfc) # expect named vector mapping KEGG id -> name

  # safe lookup: if not found, use original id or empty string
  labels <- vapply(compounds_in_graph, function(id) {
    val <- NA_character_
    if (grepl("^C", id)) {
      tmp <- compounds[compounds[[1]] == id, 2]
      val <- if (length(tmp) > 0) tmp[1] else NA_character_
    } else if (grepl("^G", id)) {
      tmp <- glycan[glycan[1] == id, 2]
      val <- if (length(tmp) > 0) tmp[1] else NA_character_
    }
    if (is.na(val)) {
      return(id)
    }

    val <- gsub(";.*", "", val) # take first name before ';'
    return(as.character(val))
  }, character(1))

  nodes_df$label[idx] <- labels
  return(nodes_df)
}
