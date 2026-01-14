#' Create a visNetwork graph from nodes and edges data frames
#' @param nodes_df Data frame of nodes.
#' @param edges_df Data frame of edges.
#' @param pathway_name Name of the pathway for the graph title.
#' @return A visNetwork object representing the graph.
#' @noRd
make_vis_graph <- function(nodes_df, edges_df, pathway_name) {
  # Different handling if no edges
  if (nrow(edges_df) == 0 || is.null(edges_df)) {
    warning("No edges in graph.")
    v <- visNetwork::visNetwork(nodes = nodes_df, main = pathway_name) # if graph has no edges
  } else {
    v <- visNetwork::visNetwork(nodes = nodes_df, edges = edges_df, main = pathway_name) # if graph has edges
  }

  v <- visNetwork::visPhysics(v, enabled = FALSE)

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
  # Map arrow.mode from igraph to visNetwork arrows
  # igraph arrow.mode: 0 = none, 1 = back, 2 = to, 3 = tee, 4 = both (etc)
  edges_df$arrows <- switch_arrow <- function(mode) {
    if (is.na(mode)) {
      return("")
    }
    switch(as.character(mode),
      "0" = "",
      "1" = "from",
      "2" = "to",
      "3" = "to", 
      "4" = "to;from",
      "" # default
    )
  }
  edges_df$arrows <- vapply(edges_df$arrow.mode, switch_arrow, character(1))

  # Map lty (igraph) to dashes (visNetwork)
  # lty = 1 solid, lty = 2 dashed
  edges_df$dashes <- ifelse(is.na(edges_df$lty), FALSE, edges_df$lty != 1)

  # Ensure color and label exist
  if (!"color" %in% names(edges_df)) edges_df$color <- "gray"
  if (!"label" %in% names(edges_df)) edges_df$label <- ""

  edges_df
}

#' Map KEGG-styled nodes to visNetwork attributes
#' @param nodes_df Data frame of nodes extracted from igraph using as_data_frame(what="vertices")
#' @return nodes_df with visNetwork-compatible styling columns: shape, borderRadius, widthConstraint, heightConstraint
#' @noRd
kegg_nodes_to_visNetwork <- function(nodes_df) {
  nodes_df$borderRadius <- NA_integer_
  # for visNetwork: set borderRadius for roundrectangle
  nodes_df$borderRadius <- ifelse(nodes_df$original_shape == "roundrectangle", 10, 0)

  # Map KEGG types to shapes
  nodes_df$shape[nodes_df$original_shape == "rectangle"] <- "box"
  nodes_df$shape[nodes_df$original_shape == "circle"] <- "dot"
  nodes_df$shape[nodes_df$original_shape == "roundrectangle"] <- "box"
  nodes_df$shape[nodes_df$original_shape == "line"] <- "ellipse"

  nodes_df$widthConstraint <- nodes_df$width # I decided to keep them interpretable in the network before and convert them to visnetwork here 
  nodes_df$heightConstraint <- nodes_df$height

  return(nodes_df)
}
