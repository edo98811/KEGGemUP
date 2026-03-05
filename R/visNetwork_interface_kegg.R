#' Create a visNetwork graph from nodes and edges data frames
#' @param vertices_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return A visNetwork object representing the graph.
#' @noRd
make_vis_graph <- function(vertices_df, edges_df, pathway_name) {
  # Different handling if no edges
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = vertices_df, main = pathway_name) # if graph has no edges
  } else {
    v <- visNetwork::visNetwork(nodes = vertices_df, edges = edges_df, main = pathway_name) # if graph has edges
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)
  v <- visNetwork::visNodes(
    v,
    shape = "dot",
    widthConstraint = FALSE
  )

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
#' @param edges_df Data frame of edges extracted from igraph using as_data_frame(what="edges")
#' @return edges_df with visNetwork-compatible styling columns: color, arrows, dashes, label
#' @noRd
igraph_edges_to_visNetwork <- function(edges_df) {

  if (is.null(edges_df) || nrow(edges_df) == 0) {
    return(edges_df)
  }
  # Map arrow.mode from igraph to visNetwork arrows
  # igraph arrow.mode: 0 = none, 1 = back, 2 = to, 3 = tee, 4 = both (etc)
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

  edges_df
}

#' Map KEGG-styled nodes to visNetwork attributes
#' @param vertices_df Data frame of nodes extracted from igraph using as_data_frame(what="vertices")
#' @param scaling_factor Numeric scaling factor for node sizes (default: 1.5)
#' @return vertices_df with visNetwork-compatible styling columns: shape, borderRadius, widthConstraint, heightConstraint
#' @noRd
kegg_nodes_to_visNetwork <- function(vertices_df, scaling_factor = 1.5) {

  if (is.null(vertices_df) || nrow(vertices_df) == 0) {
    return(vertices_df)
  }

  # Set borderRadius for roundrectangle nodes
  vertices_df$borderRadius <- ifelse(vertices_df$graphics_type == "roundrectangle", 10, 0)

  # Fix position for line nodes
  vertices_df$fixed <- ifelse(vertices_df$graphics_type == "line", TRUE, FALSE)

  # Map KEGG types to shapes
  vertices_df$shape[vertices_df$graphics_type == "rectangle"] <- "box"
  vertices_df$shape[vertices_df$graphics_type == "circle"] <- "dot"
  vertices_df$shape[vertices_df$graphics_type == "roundrectangle"] <- "box"
  vertices_df$shape[vertices_df$graphics_type == "line"] <- "text"
  vertices_df$shape[vertices_df$graphics_type == "ellipse"] <- "dot"

  vertices_df$font.size[vertices_df$graphics_type == "line"] <- 7

  vertices_df$shape[vertices_df$graphics_type == "group"] <- "dot"
  vertices_df$widthConstraint <- vertices_df$width

  vertices_df$widthConstraint <- ifelse(vertices_df$shape == "dot", NA, vertices_df$width)
  vertices_df$heightConstraint <- vertices_df$height
  vertices_df$id <- as.character(vertices_df$name)
  vertices_df$borderWidth <- 2
  vertices_df$widthConstraint[vertices_df$graphics_type == "line"] <- nchar(as.character(
    vertices_df$label[vertices_df$graphics_type == "line"] )) * 4 
  vertices_df$font.background[vertices_df$graphics_type == "line"] <- "white"
  # Set border color normally black and red on hover (except for line nodes)
  vertices_df$color <- lapply(seq_len(nrow(vertices_df)), function(i) {
    border_color <-
      if (#vertices_df$graphics_type[i] == "line" || 
      vertices_df$graphics_type[i] %in% c("group", "line")) {
        "transparent"
      } else {
        "black"
      }
    list(
      background = vertices_df$vertex.color[i],
      border = border_color,
      highlight = list(border = "red")
    )
  })

  vertices_df <- scale_dimensions(vertices_df, factor = scaling_factor)
  # vertices_df <- vertices_df[order(vertices_df$label), ]

  return(vertices_df)
}
