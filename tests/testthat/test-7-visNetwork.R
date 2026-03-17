
test_that("igraph_edges_to_visNetwork correctly maps arrows, dashes, color, label", {
  # Example edge data frame from igraph
  edges_df <- data.frame(
    from = c("A", "B", "C"),
    to = c("B", "C", "A"),
    arrow.mode = c(0, 2, NA),
    lty = c(1, 2, NA),
    stringsAsFactors = FALSE
  )

  edges_mapped <- igraph_edges_to_visNetwork(edges_df)

  # Check new columns exist
  expect_true(all(c("arrows", "dashes", "color", "label") %in% names(edges_mapped)))

  # Check arrows mapping
  expect_equal(edges_mapped$arrows, c("", "to", ""))

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
    y = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  nodes_mapped <- kegg_nodes_to_visNetwork(vertices_df, visualisation_type = "standard", scaling_factor = 1)

  # Check borderRadius exists
  expect_true("borderRadius" %in% names(nodes_mapped))
  
  # Check borderRadius
  expect_equal(nodes_mapped$borderRadius, c(0, 0, 10, 0))

  # Check shape mapping
  expect_equal(nodes_mapped$shape, c("box", "dot", "box", "text"))

  # Check that roundrectangle maps to box shape with borderRadius
  expect_equal(nodes_mapped$shape[3], "box")

  expect_equal(nodes_mapped$borderRadius[3], 10)
})


test_that("igraph_edges_to_visNetwork does not crashif edges_df empty", {
  # Example edge data frame from igraph
  edges_df <- data.frame(
    from = character(0),
    to = character(0),
    arrow.mode = integer(0),
    lty = integer(0),
    stringsAsFactors = FALSE
  )

  edges_mapped <- igraph_edges_to_visNetwork(edges_df)

  expect_equal(edges_mapped, edges_df)
})