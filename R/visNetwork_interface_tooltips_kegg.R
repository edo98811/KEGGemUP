# creating graphs in visNetwork -------------------------------------------

#' Create a visNetwork graph from nodes and edges data frames
#'
#' @details
#' This function uses a custom renderer to handle all aspects that are not
#' directly handled within visNetwork with a full R-based API
#'
#' @param vertices_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#'
#' @return A visNetwork object representing the graph.
#'
#' @importFrom htmlwidgets JS
#'
#' @noRd
make_vis_graph <- function(vertices_df,
                           edges_df,
                           pathway_name) {

# if using an external file
# renderer_dep <- htmltools::htmlDependency(
#   name = "custom-renderer",
#   version = "1.0.0",
#   src = system.file("htmlwidgets", package = "yourpackage"),
#   script = "custom_renderer.js"
# )
# custom_renderer <- htmlwidgets::JS("customRenderer")

  custom_renderer <- htmlwidgets::JS("
function(params) {

  const { ctx, x, y, style, label, state } = params;
  const { selected, hover } = state;

  const size = style?.size || 80; // this is the height for rectangle nodes
  const fillColor = style?.color || '#FFFFFF';
  // make this aware to the color, passing this directly did not really work
  const isFaded = fillColor === 'rgba(200,200,200,0.4)';

  const width = size * 2.7;
  const height = size;
  return {
    drawNode: function() {
      // only draw the rectangle here
      ctx.save();

      ctx.fillStyle = fillColor;

      // assign based on the status of being marked as faded
      ctx.strokeStyle = isFaded
        ? 'rgba(120,120,120,0.45)'
        : 'rgba(0,0,0,1)';

      ctx.lineWidth = isFaded
        ? 1
        : 2;

      ctx.beginPath();
      ctx.rect(x - width/2, y - height/2, width, height);
      ctx.fill();
      ctx.stroke();
      ctx.restore();
    },

    drawExternalLabel: function() {
      // this handles the part of drawing this *on top of the rest*
      ctx.save();
      ctx.font = '14px IBM Plex Sans'; // default could also be '14px Arial'
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';

      const text = label || '';
      const tx = x;
      const ty = y - height/2 - 1;

      const metrics = ctx.measureText(text);
      const textWidth = metrics.width;
      const textHeight = 14;

      // const padX = 8;
      // const padY = 5;

      // text
      ctx.fillStyle = isFaded
        ? 'rgba(120,120,120,0.75)'
        : '#000000';
      ctx.fillText(text, tx, ty);

      ctx.restore();
    }
  }
}
")

  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    # if graph has no edges
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = vertices_df,
                                background = .color_visnetwork_background,
                                main = list(text = pathway_name,
                                            style = .font_visnetwork_title))
  } else {
    # if graph has edges
    v <- visNetwork::visNetwork(nodes = vertices_df,
                                edges = edges_df,
                                background = .color_visnetwork_background,
                                main = list(text = pathway_name,
                                            style = .font_visnetwork_title)

    )
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)
  v <- visNetwork::visNodes(
    v,
    shape = "dot",
    widthConstraint = FALSE,
    ctxRenderer = custom_renderer
  )
  # v <- htmltools::attachDependencies(v, renderer_dep)

  v <- visNetwork::visOptions(v,
    highlightNearest = list(
      enabled = FALSE,
      # degree = 2,

      hover = FALSE
    ),
    # selectedBy = "group",
    # nodesIdSelection = TRUE
  )

  v <- visNetwork::visInteraction(v,
    dragNodes = TRUE,
    multiselect = TRUE,
    selectable = TRUE
  )
  v <- visNetwork::visEvents(v,
    selectNode = "function(nodes) {
        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});
      }",
    deselectNode = "function(nodes) {
        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});
      }"
  )

  return(v)
}

#' Map igraph-styled edges to visNetwork attributes
#'
#' @param edges_df Data frame of edges extracted from igraph using
#' `as_data_frame(x, what="edges")`
#'
#' @return edges_df with visNetwork-compatible styling columns: color, arrows,
#' dashes, label
#'
#' @noRd
kegg_edges_to_visNetwork <- function(edges_df,
                                     relationships = c("all", "reactions", "relations")) {
  relationships <- match.arg(relationships)

  if (relationships == "reactions") {
    edges_df <- edges_df[edges_df$reaction_type %in%
      c(
        "reaction_product_reversible",
        "reaction_substrate_reversible",
        "reaction_product_irreversible",
        "reaction_substrate_irreversible"
      ), ]
  } else if (relationships == "relations") {
    edges_df <- edges_df[edges_df$type == "relation", ]
  }

  if (is.null(edges_df) || nrow(edges_df) == 0) {
    return(edges_df)
  }
  # Map arrow.mode from igraph to visNetwork arrows
  # igraph arrow.mode: 0 = none, 1 = back, 2 = to, 3 = tee, 4 = both (etc)

  edges_df$arrows <- vapply(edges_df$arrow.mode, switch_arrow, character(1))

  # Map lty (igraph) to dashes (visNetwork)
  # lty = 1 solid, lty = 2 dashed
  edges_df$dashes <- ifelse(is.na(edges_df$lty), FALSE, edges_df$lty != 1)
  edges_df$dashes[edges_df$lty == 5] <- TRUE
  # Ensure color and label exist
  if (!"color" %in% names(edges_df)) edges_df$color <- "gray"
  if (!"label" %in% names(edges_df)) edges_df$label <- ""
  edges_df$from <- as.character(edges_df$from)
  edges_df$to <- as.character(edges_df$to)

  return(edges_df)
}

#' Map KEGG-styled nodes to visNetwork attributes
#'
#' @param vertices_df Data frame of nodes extracted from igraph using
#' `as_data_frame(x, what = "vertices")`
#' @param scaling_factor Numeric scaling factor for node sizes
#' @param visualization_type Character vector specifying the type of
#' visualization for nodes: "standard", "positions", or "node_name"
#'
#' @return vertices_df with visNetwork-compatible styling columns: shape,
#' borderRadius, widthConstraint, heightConstraint
#'
#' @noRd
kegg_nodes_to_visNetwork <- function(vertices_df,
                                     scaling_factor,
                                     visualization_type) {
  if (is.null(vertices_df) || nrow(vertices_df) == 0) {
    return(vertices_df)
  }

  # Set borderRadius for roundrectangle nodes
  vertices_df$borderRadius <- ifelse(vertices_df$graphics_type == "roundrectangle", 10, 0)

  # Fix position for line nodes
  vertices_df$fixed <- ifelse(vertices_df$graphics_type == "line", TRUE, FALSE)

  # Map KEGG types to shapes
  vertices_df$shape[vertices_df$graphics_type == "rectangle"] <- "custom"
  vertices_df$size[vertices_df$graphics_type == "rectangle"] <-
    vertices_df$height[vertices_df$graphics_type == "rectangle"]
  vertices_df$shape[vertices_df$graphics_type == "circle"] <- "dot"
  vertices_df$shape[vertices_df$graphics_type == "roundrectangle"] <- "box"
  vertices_df$shape[vertices_df$graphics_type == "line"] <- "text"
  vertices_df$shape[vertices_df$graphics_type == "ellipse"] <- "dot"
  # vertices_df[vertices_df$shape == "box", "margin.top"] <- vertices_df[vertices_df$shape == "box", "height"] * 0.5 + 10 # label will appear above the box

  vertices_df$font.size[vertices_df$graphics_type == "line"] <- .font_node_size

  vertices_df$shape[vertices_df$graphics_type == "group"] <- "dot"
  vertices_df$widthConstraint <- vertices_df$width

  vertices_df$widthConstraint <- ifelse(vertices_df$shape == "dot", NA, vertices_df$width)
  vertices_df$heightConstraint <- vertices_df$height
  vertices_df$id <- as.character(vertices_df$name)
  vertices_df$borderWidth <- 2
  vertices_df$widthConstraint[vertices_df$graphics_type == "line"] <- nchar(as.character(
    vertices_df$label[vertices_df$graphics_type == "line"]
  )) * 4
  vertices_df$font.background[vertices_df$graphics_type == "line"] <- .color_nodelabel_background
  # Set border color normally black and red on hover (except for line nodes)
  vertices_df$color <- lapply(seq_len(nrow(vertices_df)), function(i) {
    border_color <-
      if ( # vertices_df$graphics_type[i] == "line" ||
        vertices_df$graphics_type[i] %in% c("group", "line")) {
        .color_nodeborder_groupline
      } else {
        .color_nodeborder_default
      }
    list(
      background = vertices_df$vertex.color[i],
      border = border_color,
      highlight = list(border = .color_nodeborder_highlighted)
    )
  })

  if (visualization_type == "positions") {
    vertices_df$label <- paste0(
      "x:", round(vertices_df$x, 1), "\n",
      "y:", round(vertices_df$y, 1)
    )
  } else if (visualization_type == "node_name") {
    vertices_df$label <- vertices_df$name
  } else if (visualization_type == "node_size") {
    vertices_df$label <- paste0(
      "w:", round(vertices_df$width, 1), "\n",
      "h:", round(vertices_df$height, 1)
    )
  }

  vertices_df <- scale_dimensions(vertices_df, factor = scaling_factor)
  # vertices_df <- vertices_df[order(vertices_df$label), ]

  return(vertices_df)
}


#' Handles the arrow switch
#'
#' @param mode The mode as a character
#'
#' @returns The label for the arrows
#'
#' @noRd
switch_arrow <- function(mode) {
  if (is.na(mode)) {
    return("")
  }
  switch(as.character(mode),
         "0" = "",
         "1" = "from",
         "2" = "to",
         "3" = "to;from",
         ""
  )
}


# tooltip functions -------------------------------------------------------

#' Add tooltips to nodes for visNetwork visualization.
#'
#' @param vertices_df Data frame of nodes with columns: KEGG, label, source,
#' value.
#'
#' @return vertices_df with added 'title' column for tooltips.
#'
#' @details The tooltip includes a button to the specific KEGG entry page.
#' If multiple KEGG IDs are present, they are concatenated with '+' in the URL.
#' It also adds information about the node name, source of differential
#' expression data, and value.
#'
#' @noRd
add_node_tooltip <- function(vertices_df) {
  vertices_df$title <- paste0(
    "<div style='background:white; padding:4px;'>",
    vapply(
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
    ),
    "</div>"
  )

  return(vertices_df)
}



#' Add tooltips to edges for visNetwork visualization.
#'
#' @param edges_df Data frame of edges with columns: relation_subtype, type,
#' label.
#'
#' @return edges_df with added 'title' column for tooltips.
#'
#' @noRd
add_edge_tooltip <- function(edges_df) {
  base_url <- "https://www.genome.jp/dbget-bin/www_bget?"

  button_html_reaction <- ifelse(
    is.na(edges_df$reaction_name),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", base_url, edges_df$reaction_name, "' target='_blank'>",
      "<button type='button' style='", .style_button_edge,"'>",
      paste0(edges_df$reaction_name, " on KEGG"),
      "</button></a></div>"
    )
  )
  # Initialize the title column
  edges_df$title <- ""

  # Relation edges
  idx_relation <- edges_df$type == "relation"
  edges_df$title[idx_relation] <- paste0(
    "<div style='background:white; padding:4px;'>",
    "<h4 style='text-align: center;'>", edges_df$type[idx_relation], "</h4>",
    "<table>",
    "<tr><th align='left'>Type: </th><td>", edges_df$relation_type[idx_relation], "</td></tr>",
    "<tr><th align='left'>Subtype: </th><td>", edges_df$relation_subtype_name[idx_relation], "</td></tr>",
    # "<tr><th align='left'>Id: </th><td>", edges_df$relation_subtype_value[idx_relation], "</td></tr>",
    "</table>",
    "</div>"
  )

  # Reaction edges
  idx_reaction <- edges_df$type == "reaction"
  edges_df$title[idx_reaction] <- paste0(
    "<div style='background:white; padding:4px;'>",
    "<h4 style='text-align: center;'>", edges_df$reaction_name[idx_reaction], "</h4>",
    "<table>",
    "<tr><th align='left'>Name: </th><td>", edges_df$reaction_name[idx_reaction], "</td></tr>",
    "<tr><th align='left'>Type: </th><td>", edges_df$reaction_type[idx_reaction], "</td></tr>",
    # "<tr><th align='left'>ID: </th><td>", edges_df$reaction_id[idx_reaction], "</td></tr>",
    "</table>",
    button_html_reaction[idx_reaction],
    "</div>"
  )

  # Line edges
  idx_line <- edges_df$type == "line"
  edges_df$title[idx_line] <- paste0(
    "<div style='background:white; padding:4px;'>",
    "<h4 style='text-align: center;'>", edges_df$type[idx_line], "</h4>",
    "<table>",
    "<tr><th align='left'>Name: </th><td>", edges_df$reaction_name[idx_line], "</td></tr>",
    "</table>",
    button_html_reaction[idx_line],
    "</div>"
  )

  return(edges_df)
}


# tooltips for regular/other/group nodes ----------------------------------

#' Add tooltips to regular nodes
#'
#' @param vertices_df Data frame of nodes
#'
#' @return HTML string for the tooltip
#'
#' @noRd
regular_node_html <- function(vertices_df) {
  button_html_name <- ifelse(
    is.na(vertices_df$link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", vertices_df$link, "' target='_blank'>",
      "<button type='button' style='", .style_button_node_name,"'>",
      paste0(vertices_df$label, " on KEGG"),
      "</button></a></div>"
    )
  )

  button_html_reaction <- ifelse(
    is.na(vertices_df$reaction_link),
    "",
    paste0(
      "<div style='text-align:center; margin-top:5px;'>",
      "<a href='", vertices_df$reaction_link, "' target='_blank'>",
      "<button type='button' style='", .style_button_node_reaction, "'>",
      paste0(vertices_df$reaction_label, " on KEGG"),
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



#' Add tooltips to other nodes
#'
#' @param vertices_df Data frame of nodes
#'
#' @return HTML string for the tooltip
#'
#' @noRd
other_node_html <- function(vertices_df) {
  link_to_pathway <- paste0("https://www.kegg.jp/dbget-bin/www_bget?",
                            gsub("^path:", "", vertices_df$KEGG))

  paste0(
    "<h4 style='text-align: center;'>", vertices_df$label, "</h4>",
    "<table>",
    "<tr><th align='left'>Name </th><td>",
    ifelse(is.na(vertices_df$KEGG), "N/A", vertices_df$KEGG),
    "</td></tr>",
    # "<tr><th align='left'>ID </th><td>",
    # ifelse(is.na(vertices_df$name), "N/A", vertices_df$name),
    # "</td></tr>",
    "</table>",
    "<div style='text-align:center; margin-top:5px;'>",
    "<a href='", link_to_pathway, "' target='_blank'>",
    "<button type='button' style='", .color_button_node_pathways, "'>",
    paste0(vertices_df$KEGG, " on KEGG"),
    "</button></a></div>"
  )
}

#' Add tooltips to group nodes
#'
#' @param vertices_df Data frame of nodes
#'
#' @return HTML string for the tooltip
#'
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


