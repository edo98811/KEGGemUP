test_that("test graph is correctly generated", {
  fake_bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kegg_compounds = function(...) real_compounds,
    get_kegg_glycans = function(...) real_glycans,
    get_and_cache_kgml = function(...) kgml_path,
    parse_kgml_entries = function(...) expected_nodes,
    parse_kgml_relations = function(...) expected_edges,
    BiocFileCache = function(...) fake_bfc,
    {
      graph <- kegg_to_graph("hsa00010") # example
    }
  )


  expect_true(inherits(graph, "igraph"))
  expect_equal(igraph::vcount(graph), nrow(expected_nodes))
  expect_equal(igraph::ecount(graph), nrow(expected_edges))


  with_mocked_bindings(
    get_kegg_compounds = function(...) real_compounds,
    get_kegg_glycans = function(...) real_glycans,
    get_and_cache_kgml = function(...) kgml_path,
    parse_kgml_entries = function(...) expected_nodes,
    parse_kgml_relations = function(...) expected_edges,
    BiocFileCache = function(...) fake_bfc,
    {
      graph <- kegg_to_graph("hsa00010", return_type = "visNetwork") # example
    }
  )


  expect_true(inherits(graph, "visNetwork"))
})


test_that("results are mapped on test graph", {
  fake_bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kegg_compounds = function(...) real_compounds,
    get_kegg_glycans = function(...) real_glycans,
    get_and_cache_kgml = function(...) kgml_path,
    parse_kgml_entries = function(...) expected_nodes,
    parse_kgml_relations = function(...) expected_edges,
    BiocFileCache = function(...) fake_bfc,
    {
      graph <- kegg_to_graph("hsa00010", return_type = "igraph") # example
    }
  )

  graph_output <- suppressMessages(map_results_to_graph(graph, all_de_test_lists$genes_metabolites))
  expect_true(inherits(graph_output, "visNetwork"))

  graph_output <- suppressMessages(map_results_to_graph(graph, nodes_A_df, return_type = "visNetwork", feature_column = "KEGGID", value_column = "log2FoldChange"))
  expect_true(inherits(graph_output, "visNetwork"))
  nodes_df <- graph_output$x$nodes
  valid_shapes <- c("dot", "box")
  expect_true(all(nodes_df$shape %in% valid_shapes))

  graph_output <- suppressMessages(map_results_to_graph(graph, nodes_A_df, return_type = "igraph", feature_column = "KEGGID", value_column = "log2FoldChange"))
  expect_true(inherits(graph_output, "igraph"))
  vertex_shapes <- igraph::V(graph_output)$shape
  expect_true(all(vertex_shapes %in% c("circle", "rectangle")))
})
