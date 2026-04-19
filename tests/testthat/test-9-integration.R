test_that("workflow", {
  expect_warning(
    g_test_01 <- create_kegg_graph(
      pathway_id = "hsa00001",
      kgml_file = kgml_path_01
    ), "Could not retrieve pathway name for ID: hsa00001"
  )

  expect_s3_class(g_test_01, "igraph")

  suppressMessages(
    expect_warning(
      g_test_01_mapped <- map_results_to_graph(
        g = g_test_01,
        de_results = de_results_list
      ), "Some nodes had multiple matching IDs"
    )
  )

  expect_s3_class(g_test_01_mapped, "igraph")
  expect_true(all(c("vertex.color", "de_value", "de_name") %in% names(vertex_attr(g_test_01_mapped))))
  expect_true("legend_plots" %in% graph_attr_names(g_test_01_mapped))
  expect_type(graph_attr(g_test_01_mapped, "legend_plots"), "list")
  expect_true(all(vapply(graph_attr(g_test_01_mapped, "legend_plot"),
                         function(x) inherits(x, "gtable"),
                         FUN.VALUE = logical(1))))

  vis_graph_01 <- render_kegg_graph(g_test_01_mapped)
  expect_s3_class(vis_graph_01, "visNetwork")

  graph_01_subset <- subset_kegg_graph(g_test_01_mapped, ids_to_include = c("C00001", "C00002"))

  highlighted_graph <- highlight_kegg_graph(g_test_01_mapped, ids_to_highlight = c("C00001", "C00002"))

  expect_lt(length(V(graph_01_subset)), length(V(g_test_01_mapped)))
  expect_equal(length(V(highlighted_graph)), length(V(g_test_01_mapped)))


})
