test_that("map_results_to_graph correctly maps all DE test lists onto igraph", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    # Use the preloaded reference graph
    g <- reference_graph

    # Map results
    mapped_graph <- map_results_to_graph(
      g,
      de_list,
      feature_column = "KEGG_ids",
      value_column = "log2FoldChange",
      palette = "RdBu",
      verbose = FALSE
    )

    expect_s3_class(mapped_graph, "igraph", info = paste0(test_name, ": not an igraph"))

    # Check key vertex attributes
    for (attr in c("de_value", "de_source", "color")) {
      expect_true(attr %in% igraph::vertex_attr_names(mapped_graph),
        info = paste0(test_name, ": missing vertex attribute ", attr)
      )
    }

    # Values should be numeric or NA
    de_vals <- igraph::vertex_attr(mapped_graph, "de_value")
    expect_true(all(is.numeric(de_vals) | is.na(de_vals)),
      info = paste0(test_name, ": de_value contains non-numeric values")
    )

    # Color should be hex codes or NA
    colors <- igraph::vertex_attr(mapped_graph, "color")
    non_na_colors <- colors[!is.na(colors)]
    expect_true(all(grepl("^#([A-Fa-f0-9]{6})$", non_na_colors)),
      info = paste0(test_name, ": invalid hex color codes detected")
    )
  })
})

test_that("kegg_to_graph returns correct igraph using local KGML", {
  # Run function
  g <- kegg_to_graph(
    pathway_id = "hsa00010",
    scaling_factor = 1.0,
    simplified_graph = TRUE,
    kgml_file = kgml_path_01
  )

  # Basic checks
  expect_s3_class(g, "igraph")

  # Compare structure: number of nodes and edges
  expect_equal(vcount(g), vcount(ref_graph))
  expect_equal(ecount(g), ecount(ref_graph))

  # Compare vertex names
  expect_setequal(V(g)$name, V(ref_graph)$name)

  expect_equal(g, expected_graphs$g_test_01)

  # Run function
  g <- kegg_to_graph(
    pathway_id = "hsa00010",
    scaling_factor = 1.0,
    simplified_graph = TRUE,
    kgml_file = kgml_path_02
  )

  # Basic checks
  expect_s3_class(g, "igraph")

  # Compare structure: number of nodes and edges
  expect_equal(vcount(g), vcount(ref_graph))
  expect_equal(ecount(g), ecount(ref_graph))

  # Compare vertex names
  expect_setequal(V(g)$name, V(ref_graph)$name)
  expect_equal(g, expected_graphs$g_test_02)

  # Check ids_for_mapping attribute
  ids_mapping <- igraph::vertex_attr(g, "ids_for_mapping")
  expect_true(is.character(ids_mapping))
  expect_true(!all(nchar(ids_mapping[!is.na(ids_mapping)]) == 0))
})


test_that("plot_visNetwork_kegg works correctly", {
  # Your input igraph
  g_input <- reference_graph

  # Your expected visNetwork object
  v_expected <- reference_visnetwork

  # Run the function
  v_actual <- plot_visNetwork_kegg(g_input)

  # Basic type check
  expect_s3_class(v_actual, "visNetwork")

  # Check that key components exist
  expect_true(all(c("nodes", "edges") %in% names(v_actual)))

  # Compare nodes data frames
  expect_equal(
    v_actual$nodes[order(v_actual$nodes$id), ],
    v_expected$nodes[order(v_expected$nodes$id), ],
    tolerance = 1e-8
  )

  # Compare edges data frames
  expect_equal(
    v_actual$edges[order(v_actual$edges$from, v_actual$edges$to), ],
    v_expected$edges[order(v_expected$edges$from, v_expected$edges$to), ],
    tolerance = 1e-8
  )

  # Optionally compare graph-level attributes if present
  if ("graph_attr" %in% names(v_actual) && "graph_attr" %in% names(v_expected)) {
    expect_equal(v_actual$graph_attr, v_expected$graph_attr)
  }
})
