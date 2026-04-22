
test_that("kegg_edges_to_visNetwork correctly maps arrows, dashes, color, label", {
  # Example edge data frame from igraph
  edges_df <- data.frame(
    from = c("A", "B", "C"),
    to = c("B", "C", "A"),
    arrow.mode = c(0, 2, 1),
    lty = c(1, 2, NA)
  )

  edges_mapped <- kegg_edges_to_visNetwork(edges_df)

  # Check new columns exist
  expect_true(all(c("arrows", "dashes", "color", "label") %in% names(edges_mapped)))

  # Check arrows mapping
  expect_equal(edges_mapped$arrows, c("", "to", "from"))

  # Check dashes mapping
  expect_equal(edges_mapped$dashes, c(FALSE, TRUE, FALSE))

  # Check default color and label
  expect_equal(edges_mapped$color, rep("gray", 3))
  expect_equal(edges_mapped$label, rep("", 3))
})

test_that("kegg_nodes_to_visNetwork correctly maps shapes and borderRadius", {
  # Example nodes data frame
  vertices_df <- data.frame(
    name = c("n1", "n2", "n3", "n4"),
    graphics_type = c("rectangle", "circle", "roundrectangle", "line"),
    width = c(50, 40, 30, 20),
    height = c(20, 30, 40, 50),
    shape = NA_character_,
    label = c("Node1", "Node2", "Node3", "Node4"),
    x = c(1, 2, 3, 4),
    y = c(1, 2, 3, 4)
  )

  nodes_mapped <- kegg_nodes_to_visNetwork(vertices_df, visualisation_type = "standard", scaling_factor = 1)

  # Check borderRadius exists
  expect_true("borderRadius" %in% names(nodes_mapped))

  # Check borderRadius
  expect_equal(nodes_mapped$borderRadius, c(0, 0, 10, 0))

  # Check shape mapping
  expect_equal(nodes_mapped$shape, c("custom", "dot", "box", "text"))

  # Check that roundrectangle maps to box shape with borderRadius
  expect_equal(nodes_mapped$shape[3], "box")
  expect_equal(nodes_mapped$borderRadius[3], 10)
})


test_that("kegg_edges_to_visNetwork does not crash if edges_df empty", {
  # Example edge data frame from igraph
  edges_df <- data.frame(
    from = character(0),
    to = character(0),
    arrow.mode = integer(0),
    lty = integer(0)
  )

  edges_mapped <- kegg_edges_to_visNetwork(edges_df)

  expect_equal(edges_mapped, edges_df)
})

test_that("test_edge_tooltip works correctly", {
  g_input <- expected_output$g_test_01_mapped
  edges_df <- as_data_frame(g_input, what = "edges")
  edges_df <- add_edge_tooltip(edges_df)

  rel_titles <- edges_df[edges_df$type == "relation", "title", drop = FALSE]
  expect_true(all(grepl("relation", rel_titles$title, ignore.case = TRUE)))

  react_titles <- edges_df[edges_df$type == "reaction", "title", drop = FALSE]
  expect_true(all(grepl("reaction", react_titles$title, ignore.case = TRUE)))
  expect_true(all(grepl("<button", react_titles$title, ignore.case = TRUE)))

  line_titles <- edges_df[edges_df$type == "line", "title", drop = FALSE]
  expect_true(all(grepl("line", line_titles$title, ignore.case = TRUE)))
})

test_that("test_node_tooltip works correctly", {

  g_input <- expected_output$g_test_01_mapped
  nodes_df <- as_data_frame(g_input, what = "vertices")
  nodes_df <- add_node_tooltip(nodes_df)

  nodes_df_group <- nodes_df[nodes_df$type == "group", "title", drop = FALSE]
  nodes_df_compound <- nodes_df[nodes_df$type == "compound", "title", drop = FALSE]
  nodes_df_ortholog <- nodes_df[nodes_df$type == "ortholog", "title", drop = FALSE]
  nodes_df_gene <- nodes_df[nodes_df$type == "gene", "title", drop = FALSE]
  nodes_df_map <- nodes_df[nodes_df$type == "map", "title", drop = FALSE]

  expect_true(all(grepl("group", nodes_df_group$title, ignore.case = TRUE)))
  expect_true(all(grepl("KEGG ID", nodes_df_compound$title, ignore.case = TRUE)))
  expect_true(all(grepl("KEGG ID", nodes_df_ortholog$title, ignore.case = TRUE)))
  expect_true(all(grepl("KEGG ID", nodes_df_gene$title, ignore.case = TRUE)))
  expect_true(all(grepl("map", nodes_df_map$title, ignore.case = TRUE)))

  expect_true(all(grepl("<button", nodes_df_compound$title, ignore.case = TRUE)))
  expect_true(all(grepl("<button", nodes_df_ortholog$title, ignore.case = TRUE)))
  expect_true(all(grepl("<button", nodes_df_gene$title, ignore.case = TRUE)))
  # expect_true(all(grepl("<button", nodes_df_map$title, ignore.case = TRUE))) to add
})
