# test_that("test graph is correctly generated", {
#   # Load mapping and mock get_compounds to return it
#   # 4. Fake BFC
#   fake_bfc <- BiocFileCache(tempfile(), ask = FALSE)

#   # 5. Mock functions with correct return TYPES
#   with_mocked_bindings(
#     get_kegg_compounds = mock(real_compounds),
#     get_kegg_glycans = mock(real_glycans),
#     download_kgml = mock(kgml_path),
#     parse_kgml_entries = mock(expected_nodes),
#     parse_kgml_relations = mock(expected_edges),
#     BiocFileCache = mock(fake_bfc),
#     {
#       g <- kegg_to_graph("hsa00010") # example
#     }
#   )

#   # 6. Assertions
#   expect_true(inherits(g_expected, "igraph"))
#   expect_equal(igraph::vcount(g), nrow(expected_nodes))
#   expect_equal(igraph::ecount(g), nrow(expected_edges))

#   # 5. Mock functions with correct return TYPES
#   with_mocked_bindings(
#     get_kegg_compounds = mock(real_compounds),
#     get_kegg_glycans = mock(real_glycans),
#     download_kgml = mock(kgml_path),
#     parse_kgml_entries = mock(expected_nodes),
#     parse_kgml_relations = mock(expected_edges),
#     BiocFileCache = mock(fake_bfc),
#     {
#       g <- kegg_to_graph("hsa00010", return_type = "visNetwork") # example
#     }
#   )

#   # 6. Assertions
#   expect_true(inherits(g_expected, "visNetwork"))
# })
# test_that("values are correctly mapped to test graph", {

# })
