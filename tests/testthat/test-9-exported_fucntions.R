# test_that("map_results_to_graph correctly maps DE results onto igraph", {
#   # Use the preloaded reference graph
#   g <- expected_graphs$g_test_02

#   # Map results
#   expect_warning(
#     mapped_graph <-
#       map_results_to_graph(
#         g,
#         de_results_list,
#         feature_column = "KEGG_ids",
#         value_column = "log2FoldChange",
#         palette = "RdBu"
#       ),
#     "Some nodes had multiple matching IDs; only the first match was used for de_value/de_source."
#   )

#   expect_s3_class(mapped_graph, "igraph")

#   # Check key vertex attributes
#   for (attr in c("de_value", "de_source", "vertex.color")) {
#     expect_true(
#       attr %in% igraph::vertex_attr_names(mapped_graph),
#       info = paste0("missing vertex attribute ", attr)
#     )
#   }

#   # Values should be numeric or NA
#   de_vals <- igraph::vertex_attr(mapped_graph, "de_value")
#   expect_true(
#     all(is.numeric(de_vals) | is.na(de_vals)),
#     info = "de_value contains non-numeric values"
#   )

#   # vertex.color should be non NULL
#   vertex.colors <- igraph::vertex_attr(mapped_graph, "vertex.color")

#   # test for correctness
#   expect_true(any(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", vertex.colors)))
#   expect_true(any(vertex.colors != "white"))
# })


# test_that("kegg_to_graph returns correct igraph using local KGML", {
#   # Run function
#   expect_warning(
#     g <-
#       kegg_to_graph(
#         pathway_id = "hsa00001",
#         kgml_file = kgml_path_01
#       ),
#     "Could not retrieve pathway name for ID: hsa00001"
#   )

#   expect_s3_class(g, "igraph")

#   # Run function
#   expect_warning(
#     g <-
#       kegg_to_graph(
#         pathway_id = "hsa00001",
#         kgml_file = kgml_path_02
#       ),
#     "Could not retrieve pathway name for ID: hsa00001"
#   )

#   # Basic checks
#   expect_s3_class(g, "igraph")

#   # Check ids_for_mapping attribute
#   ids_mapping <- igraph::vertex_attr(g, "ids_for_mapping")
#   expect_true(is.character(ids_mapping))
#   expect_true(!all(nchar(ids_mapping[!is.na(ids_mapping)]) == 0))
# })

# test_that("make_kegg_visNetwork works correctly", {
#   # Your input igraph
#   g_input <- expected_graphs$g_test_01_mapped

#   # Your expected visNetwork object
#   v_expected <- expected_graphs$visNetwork_test_01_mapped

#   # Run the function
#   v_actual <- make_kegg_visNetwork(g_input)

#   # Basic type check
#   expect_s3_class(v_actual, "visNetwork")

#   # Compare to expected visNetwork object
#   expect_equal(v_actual$x$nodes, v_expected$x$nodes)
#   expect_equal(v_actual$x$edges, v_expected$x$edges)

#   # Your input igraph
#   g_input <- expected_graphs$g_test_02_mapped

#   # Your expected visNetwork object
#   v_expected <- expected_graphs$visNetwork_test_02_mapped

#   # Run the function
#   v_actual <- make_kegg_visNetwork(g_input)

#   # Basic type check
#   expect_s3_class(v_actual, "visNetwork")

#   # Compare to expected visNetwork object
#   expect_equal(v_actual$x$nodes, v_expected$x$nodes)
#   expect_equal(v_actual$x$edges, v_expected$x$edges)

#   # Compare graph-level attributes if present
#   if ("graph_attr" %in% names(v_actual) && "graph_attr" %in% names(v_expected)) {
#     expect_equal(v_actual$graph_attr, v_expected$graph_attr)
#   }
# })
