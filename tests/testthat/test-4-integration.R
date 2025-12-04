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

  graph_output <- map_results_to_nodes(graph, all_de_test_lists$genes_metabolites)
  expect_true(inherits(graph_output, "visNetwork"))

  graph_output <- map_results_to_nodes(graph, nodes_A_df, return_type = "visNetwork", feature_column = "KEGGID", value_column = "log2FoldChange")
  expect_true(inherits(graph_output, "visNetwork"))
})
