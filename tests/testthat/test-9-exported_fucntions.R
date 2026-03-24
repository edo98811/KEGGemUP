test_that("map_results_to_graph correctly maps DE results onto igraph", {
  # Use the preloaded reference graph
  g <- expected_graphs$g_test_01

  # Map results
  suppressMessages(
    expect_warning(
      mapped_graph <-
        map_results_to_graph(
          g,
          de_results_list
        ),
      "Some nodes had multiple matching IDs; only the first match was used for de_value/de_source."
    )
  )
  expect_s3_class(mapped_graph, "igraph")

  # Check key vertex attributes
  for (attr in c("de_value", "de_source", "vertex.color")) {
    expect_true(
      attr %in% igraph::vertex_attr_names(mapped_graph),
      info = paste0("missing vertex attribute ", attr)
    )
  }

  # Values should be numeric or NA
  de_vals <- igraph::vertex_attr(mapped_graph, "de_value")
  expect_true(
    all(is.numeric(de_vals) | is.na(de_vals)),
    info = "de_value contains non-numeric values"
  )

  # vertex.color should be non NULL
  vertex.colors <- igraph::vertex_attr(mapped_graph, "vertex.color")

  # test for correctness of colors
  expect_true(any(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", vertex.colors)))
  expect_true(any(vertex.colors != "white"))
})

test_that("map_results_to_graph does not accept a non igraph object", {
  expect_error(
    map_results_to_graph(
      g = data.frame(name = c("A", "B"), label = c("A", "B")),
      de_results_list
    ),
    "Input graph 'g' must be an igraph object."
  )
})

test_that("kegg_to_graph returns correct igraph using local KGML", {
  # Run function
  expect_warning(
    g <-
      kegg_to_graph(
        pathway_id = "hsa00001",
        kgml_file = kgml_path_01
      ),
    "Could not retrieve pathway name for ID: hsa00001"
  )

  expect_s3_class(g, "igraph")

  # Check ids_for_mapping attribute
  ids_mapping <- igraph::vertex_attr(g, "ids_for_mapping")
  expect_true(is.character(ids_mapping))
  expect_true(!all(nchar(ids_mapping[!is.na(ids_mapping)]) == 0))
})

test_that("kegg_to_graph handles missing KGML file", {
  expect_warning(
    expect_error(
      kegg_to_graph(pathway_id = "hsa00001", kgml_file = "non_existent_file.xml"),
      regexp = "Failed to read XML file"
    )
  )
})

test_that("make_kegg_visNetwork works correctly", {
  # Your input igraph
  g_input <- expected_graphs$g_test_01_mapped

  # Your expected visNetwork object
  v_expected <- expected_graphs$visNetwork_test_01_mapped

  # Run the function
  v_actual <- make_kegg_visNetwork(g_input)

  # Basic type check
  expect_s3_class(v_actual, "visNetwork")

  # Compare to expected visNetwork object
  expect_equal(v_actual$x$nodes, v_expected$x$nodes)
  expect_equal(v_actual$x$edges, v_expected$x$edges)
})

test_that("make_graph_subset works correctly", {
  ids_to_subset <- c("C00001", "C00002")

  g_input <- expected_graphs$g_test_01_mapped
  g_subset <- make_graph_subset(g_input, ids_to_subset)

  # Check that the subset graph only contains the specified nodes
  subset_node_names <- igraph::vertex_attr(g_subset, "ids_for_mapping")
  expect_true(all(ids_to_subset %in% subset_node_names))

  # Check that the edges in the subset graph only connect the specified nodes
  edge_list <- igraph::as_data_frame(g_subset, what = "edges")
  vertex_names <- igraph::vertex_attr(g_subset, "name")
  expect_true(all(edge_list$from %in% vertex_names))
  expect_true(all(edge_list$to %in% vertex_names))

})

test_that("make_graph_subset throws error for non-igraph input", {
  expect_error(
    make_graph_subset(g = data.frame(name = c("A", "B")), ids_to_include = c("A")),
    "Input graph 'g' must be an igraph object."
  )
})

test_that("make_graph_subset throws error for empty or wrong ids_to_include", {
  g_input <- expected_graphs$g_test_01_mapped
  expect_error(
    make_graph_subset(g = g_input, ids_to_include = character(0)),
    "ids_to_include must be a non-empty character vector of KEGG IDs."
  )

})

test_that("make_graph_subset returns original graph if no matching nodes found", {
  g_input <- expected_graphs$g_test_01_mapped
  expect_warning(
    make_graph_subset(g = g_input, ids_to_include = c("NON_EXISTENT_ID")),
    "No matching nodes found for the provided KEGG IDs."
  )
})

test_that("highlight_graph_subset works correctly", {
  ids_to_highlight <- c("C00001", "C00002")

  g_input <- expected_graphs$g_test_01_mapped
  g_highlighted <- highlight_graph_subset(g_input, ids_to_highlight)

  vertices_df <- igraph::as_data_frame(g_highlighted, what = "vertices")
  mapping <- make_mapping_df(vertices_df)
  nodes_to_highlight <- mapping$matched_id %in% ids_to_highlight
  # Check that the highlighted graph is still an igraph
  nodes_to_highlight <- unique(
    mapping[mapping$matched_id %in% ids_to_highlight, "name"]
  )
  nodes_to_highlight <- nodes_to_highlight[!is.na(nodes_to_highlight)]

  # Check that highlighted nodes have the correct color
  vertex_colors <- igraph::vertex_attr(g_highlighted, "vertex.color")
  highlighted_node_indices <- which(igraph::vertex_attr(g_highlighted, "name") %in% nodes_to_highlight)
  non_highlighted_node_indices <- setdiff(seq_along(vertex_colors), highlighted_node_indices)
  expect_true(all(vertex_colors[non_highlighted_node_indices] == "rgba(200,200,200,0.4)"))
})

test_that("highlight_graph_subset throws error for non-igraph input", {
  expect_error(
    highlight_graph_subset(g = data.frame(name = c("A", "B")), ids_to_highlight = c("A")),
    "Input graph 'g' must be an igraph object."
  )
})

test_that("highlight_graph_subset throws error for empty or wrong ids_to_highlight", {
  g_input <- expected_graphs$g_test_01_mapped
  expect_error(
    highlight_graph_subset(g = g_input, ids_to_highlight = character(0)),
    "ids_to_highlight must be a non-empty character vector of KEGG IDs."
  )
})

test_that("highlight_graph_subset returns original graph if no matching nodes found", {
  g_input <- expected_graphs$g_test_01_mapped
  expect_warning(
    highlight_graph_subset(g = g_input, ids_to_highlight = c("NON_EXISTENT_ID")),
    "No matching nodes found for the provided KEGG IDs."
  )
})