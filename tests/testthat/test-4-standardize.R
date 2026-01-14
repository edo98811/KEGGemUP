# node/edge mapping
test_standardize_network_basic <- function() {
  # Create example graph
  g <- graph_from_data_frame(
    data.frame(from = c(1, 2), to = c(2, 3), weight = c(0.5, 1.2), kegg_edge = c("e1","e2")),
    vertices = data.frame(
      name = 1:3,
      KEGG = c("K00001", "K00002", "K00003"),
      kegg_node_val = c(10, 20, 30)
    ),
    directed = TRUE
  )

  g_std <- standardize_network(g, kegg_to_general_edge_map(), kegg_to_general_node_map(), node_defaults(), edge_defaults())

  nodes_df <- as_data_frame(g_std, what = "vertices")
  edges_df <- as_data_frame(g_std, what = "edges")

  # Check node mapping
  stopifnot(all(names(node_defaults()) %in% names(nodes_df)))
  stopifnot(nodes_df$value[1] == 10 && nodes_df$type[1] == "gene")

  # Check edge mapping
  stopifnot(all(names(edge_defaults()) %in% names(edges_df)))
  stopifnot(edges_df$edge_id[1] == "e1")
  stopifnot(edges_df$weight[1] == 0.5)
}
