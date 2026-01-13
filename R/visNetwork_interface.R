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

convert_to_visnetwork_dfs <- function(nodes_df, edges_df) {
  
  # Shapes conversion for visNetwork
  nodes_df$shape[nodes_df$shape == "vrectangle"] <- "box"
  nodes_df$shape[nodes_df$shape == "circle"] <- "dot"

  return(list(nodes = nodes_df, edges = edges_df))
}