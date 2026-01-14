# Test 1: basic node/edge mapping
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

  node_map <- list(kegg_node_val = "value")
  edge_map <- list(kegg_edge = "edge_id")
  node_default <- list(value = NA, type = "gene")
  edge_default <- list(edge_id = NA, weight = 1)

  g_std <- standardize_network(g, node_map, edge_map, node_default, edge_default)

  nodes_df <- as_data_frame(g_std, what = "vertices")
  edges_df <- as_data_frame(g_std, what = "edges")

  # Check node mapping
  stopifnot(all(c("value", "type") %in% names(nodes_df)))
  stopifnot(nodes_df$value[1] == 10 && nodes_df$type[1] == "gene")

  # Check edge mapping
  stopifnot(all(c("edge_id", "weight") %in% names(edges_df)))
  stopifnot(edges_df$edge_id[1] == "e1")
  stopifnot(edges_df$weight[1] == 0.5)
  print("Test 1 passed")
}

# Test 2: missing KEGG columns use default values
test_standardize_network_defaults <- function() {
  g <- graph_from_data_frame(
    data.frame(from = 1, to = 2),
    vertices = data.frame(name = 1:2),
    directed = FALSE
  )

  node_map <- list(nonexistent_node = "value")
  edge_map <- list(nonexistent_edge = "edge_id")
  node_default <- list(value = 100, type = "protein")
  edge_default <- list(edge_id = "none", weight = 0)

  g_std <- standardize_network(g, node_map, edge_map, node_default, edge_default)

  nodes_df <- as_data_frame(g_std, what = "vertices")
  edges_df <- as_data_frame(g_std, what = "edges")

  stopifnot(all(nodes_df$value == 100))
  stopifnot(all(nodes_df$type == "protein"))
  stopifnot(all(edges_df$edge_id == "none"))
  stopifnot(all(edges_df$weight == 0))
  print("Test 2 passed")
}

# Test 3: simplified graph removes KEGG columns
test_standardize_network_simplified <- function() {
  g <- graph_from_data_frame(
    data.frame(from = 1, to = 2, KEGG_edge = "e1"),
    vertices = data.frame(name = 1:2, KEGG_node = "k1"),
    directed = TRUE
  )

  node_map <- list()
  edge_map <- list()
  node_default <- list(value = 1)
  edge_default <- list(weight = 0)

  g_std <- standardize_network(g, node_map, edge_map, node_default, edge_default, simplified_graph = TRUE)

  nodes_df <- as_data_frame(g_std, what = "vertices")
  edges_df <- as_data_frame(g_std, what = "edges")

  stopifnot(!"KEGG_node" %in% names(nodes_df))
  stopifnot(!"KEGG_edge" %in% names(edges_df))
  stopifnot("value" %in% names(nodes_df) && "weight" %in% names(edges_df))
  print("Test 3 passed")
}

# Test 4: graph with no edges
test_standardize_network_no_edges <- function() {
  g <- make_empty_graph(n = 2, directed = FALSE)
  node_map <- list()
  edge_map <- list()
  node_default <- list(value = 0)
  edge_default <- list(weight = 0)

  g_std <- standardize_network(g, node_map, edge_map, node_default, edge_default)

  nodes_df <- as_data_frame(g_std, what = "vertices")
  edges_df <- as_data_frame(g_std, what = "edges")

  stopifnot(nrow(nodes_df) == 2)
  stopifnot(nrow(edges_df) == 0)
  stopifnot(all(names(edges_df) == names(edge_default)))
  print("Test 4 passed")
}

# Run all tests
test_standardize_network_basic()
test_standardize_network_defaults()
test_standardize_network_simplified()
test_standardize_network_no_edges()
