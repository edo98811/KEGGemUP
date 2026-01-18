# test_that("standardize_network maps nodes and edges correctly", {
#   # Build and standardize the graph
#   g <- build_kegg_graph(kgml_path_02, "Test", bfc = bfc)
#   g_standardized <- standardize_network(
#     g,
#     kegg_to_general_node_map(),
#     kegg_to_general_edge_map(),
#     node_defaults(),
#     edge_defaults()
#   )

#   # Extract standardized nodes and edges
#   nodes_std <- as_data_frame(g_standardized, what = "vertices")
#   edges_std <- as_data_frame(g_standardized, what = "edges")

#   # Check column names
#   expected_node_cols <- c(
#     "name", "link", "ids_for_mapping", "label", "x", "y", "width", "height",
#     "group", "feature_id_1", "feature_id_2", "graphics_type", "title",
#     "de_value", "de_source", "color", "size", "fixed", "borderRadius",
#     "shape", "text"
#   )
#   expect_true(all(expected_node_cols %in% colnames(nodes_std)))

#   # Check that certain important nodes exist
#   expect_true(all(c("1", "2", "3") %in% nodes_std$name))

#   # Check specific values for a few key nodes
#   node_1 <- nodes_std[nodes_std$name == "1", ]
#   expect_equal(node_1$link, "https://www.kegg.jp/dbget-bin/www_bget?tst:1111")
#   expect_equal(node_1$label, "1111")
#   expect_equal(node_1$x, 100)
#   expect_equal(node_1$y, 100)

#   node_2 <- nodes_std[nodes_std$name == "2", ]
#   expect_equal(node_2$ids_for_mapping, "2222;3333")
#   expect_equal(node_2$label, "2222;3333")

#   # Check column names
#   expected_edge_cols <- c(
#     "from", "to", "type", "name_3", "label", "name_1", "name_2",
#     "id", "link", "directed", "color", "width", "value", "lty",
#     "arrows", "dashes", "title"
#   )
#   expect_true(all(expected_edge_cols %in% colnames(edges_std)))

#   # Check that some key edges exist
#   expect_true(any(edges_std$from == "1" & edges_std$to == "2" & edges_std$type == "PPrel"))
#   expect_true(any(edges_std$from == "2" & edges_std$to == "3" & edges_std$type == "ECrel"))

#   # Check that directed edges are TRUE
#   expect_true(all(edges_std$directed[!is.na(edges_std$directed)]))

#   # Optional: check edge widths and colors
#   expect_true(all(edges_std$color == "gray"))
#   expect_true(all(edges_std$width == 1))
# })

test_that("build_kegg_graph handles empty KGML files", {
  expect_error(
    build_kegg_graph(empty_kgml_path, "Empty_Pathway", bfc = bfc),
    regexp = "No nodes found in KGML file"
  )
})

test_that("build_kegg_graph works correctly", {
  g <- build_kegg_graph(kgml_path_01, "Test_Pathway", bfc = bfc)

  expect_true(inherits(g, "igraph"))
  nodes_df <- as_data_frame(g, what = "vertices")
  edges_df <- as_data_frame(g, what = "edges")

  expect_true(all(edges_df$from %in% nodes_df$name))
  expect_true(all(edges_df$to %in% nodes_df$name))

  expect_equal(igraph::graph_attr(g, "title"), "Test_Pathway")
  expect_equal(igraph::graph_attr(g, "type"), "KEGG_Pathway")

  expect_true(all(is.character(nodes_df$ids_for_mapping)))
  expect_true(all(nodes_df$ids_for_mapping == "" | nzchar(nodes_df$ids_for_mapping)))
})

test_that("make_igraph_graph handles edge cases", {
  nodes_df <- data.frame(name = c("1", "2"), label = c("A", "B"))
  edges_df <- data.frame(from = character(), to = character())

  expect_warning(
    g <- make_igraph_graph(nodes_df, edges_df, "Test_Pathway"),
    regexp = "No edges in graph"
  )

  expect_true(inherits(g, "igraph"))
  expect_equal(igraph::vcount(g), 2)
  expect_equal(igraph::ecount(g), 0)

  nodes_df <- data.frame(name = c("1", "2", "3"), label = c("A", "B", "C"))
  edges_df <- data.frame(from = c("1", "2"), to = c("2", "3"))

  g <- make_igraph_graph(nodes_df, edges_df, "Test_Pathway")

  expect_true(inherits(g, "igraph"))
  expect_equal(igraph::vcount(g), 3)
  expect_equal(igraph::ecount(g), 2)

  nodes_df <- data.frame(name = c("1", "2"), label = c("A", "B"))

  expect_warning(
    g <- make_igraph_graph(nodes_df, NULL, "Test_Pathway"),
    regexp = "No edges in graph"
  )

  expect_true(inherits(g, "igraph"))
  expect_equal(igraph::vcount(g), 2)
  expect_equal(igraph::ecount(g), 0)

  nodes_df <- data.frame(name = c("3", "1", "2"), label = c("Z", "A", "B"))
  edges_df <- data.frame(from = c("3", "1"), to = c("1", "2"))

  g <- make_igraph_graph(nodes_df, edges_df, "Test_Pathway")

  vertex_labels <- igraph::V(g)$label
  expect_equal(vertex_labels, sort(vertex_labels))
})


# node/edge mapping
# test_standardize_network_basic <- function() {
#   # Create example graph
#   g <- graph_from_data_frame(
#     data.frame(from = c(1, 2), to = c(2, 3), weight = c(0.5, 1.2), kegg_edge = c("e1", "e2")),
#     vertices = data.frame(
#       name = 1:3,
#       KEGG = c("K00001", "K00002", "K00003"),
#       kegg_node_val = c(10, 20, 30)
#     ),
#     directed = TRUE
#   )

#   g_std <- standardize_network(
#     g,
#     kegg_to_general_edge_map(),
#     kegg_to_general_node_map(),
#     node_defaults(),
#     edge_defaults()
#   )

#   nodes_df <- as_data_frame(g_std, what = "vertices")
#   edges_df <- as_data_frame(g_std, what = "edges")

#   # Check node mapping
#   stopifnot(all(names(node_defaults) %in% names(nodes_df)))
#   stopifnot(nodes_df$value[1] == 10 && nodes_df$type[1] == "gene")

#   # Check edge mapping
#   stopifnot(all(names(edge_defaults) %in% names(edges_df)))
#   stopifnot(edges_df$edge_id[1] == "e1")
#   stopifnot(edges_df$weight[1] == 0.5)
# }
