test_that("style_igraph_graph works correctly on reactions", {
  
  g_test_01_edges_df <- style_igraph_graph(kgml_steps$g_test_01)
  
  edges_df <- as_data_frame(g_test_01_edges_df, what = "edges")
  edges_df <- edges_df[
    edges_df$type == "reaction",
    c("arrow.mode","dashes","reaction_type","color","label")
  ]
  
  # split reaction types
  irreversible <- edges_df[grepl("_irreversible$", edges_df$reaction_type), ]
  reversible   <- edges_df[grepl("_reversible$", edges_df$reaction_type), ]
  
  expect_true(all(irreversible$color == "#FF4500"))
  expect_true(all(reversible$color == "#008000"))
  expect_true(all(irreversible$label == "→"))
  expect_true(all(reversible$label == "⇄"))
  expect_true(all(irreversible$arrow.mode %in% c(0,2)))
  expect_true(all(reversible$arrow.mode %in% c(1,2)))
  
  # arrows should appear equally often
  expect_equal(
    sum(reversible$arrow.mode == 1),
    sum(reversible$arrow.mode == 2)
  )
  # in irreversible one of the reactions does not have middle gene
  expect_false(
    sum(irreversible$arrow.mode == 1) == sum(irreversible$arrow.mode == 2)
  )
  
})

test_that("style_igraph_graph works correctly on relations", {

  g_test_01_edges_df <- style_igraph_graph(kgml_steps$g_test_01)
  edges_df <- as_data_frame(g_test_01_edges_df, what = "edges")
  edges_df <- edges_df[
    edges_df$type == "relation",
    c("arrow.mode","dashes","relation_type","color","label")
  ]
  
  # split based on types
  maplink <- edges_df[edges_df$relation_type == "maplink", ]
  others  <- edges_df[edges_df$relation_type != "maplink", ]

  expect_true(all(maplink$arrow.mode == 1))
  expect_true(all(others$arrow.mode == 2))
  expect_true(all(edges_df$color[edges_df$relation_type == "PPrel"] == "#FF6347"))
  expect_true(all(edges_df$color[edges_df$relation_type == "ECrel"] == "#8A2BE2"))
  expect_true(all(edges_df$color[edges_df$relation_type == "maplink"] == "#FF4500"))
  expect_true(all(edges_df$label[edges_df$relation_type == "PPrel"] == "PP"))
  expect_true(all(edges_df$label[edges_df$relation_type == "ECrel"] == "EC"))
  expect_true(all(edges_df$label[edges_df$relation_type == "maplink"] == "maplink"))
})


test_that("style_igraph_graph works correctly on vertices", {
  
  g_test_01_nodes_df <- style_igraph_graph(kgml_steps$g_test_01)
  nodes_df <- as_data_frame(g_test_01_nodes_df, what = "vertices")

  # Testing the properties...
  expect_equal(unique(nodes_df$shape[nodes_df$graphics_type == "rectangle"]), "vrectangle")
  expect_equal(unique(nodes_df$shape[nodes_df$graphics_type == "circle"]), "circle")
  expect_equal(unique(nodes_df$shape[nodes_df$graphics_type == "roundrectangle"]), "vrectangle")
  expect_equal(unique(nodes_df$shape[nodes_df$graphics_type == "line"]), "circle")
  expect_equal(unique(nodes_df$shape[nodes_df$graphics_type == "group"]), "circle")
  
  expect_equal(unique(nodes_df$size[nodes_df$graphics_type == "rectangle"]), 25)
  expect_equal(unique(nodes_df$size[nodes_df$graphics_type == "roundrectangle"]), 25)
  expect_equal(unique(nodes_df$size[nodes_df$graphics_type == "circle"]), 10)
  expect_equal(unique(nodes_df$size[nodes_df$graphics_type == "line"]), 1)
  expect_equal(unique(nodes_df$size[nodes_df$type == "group"]), 2)
  
  expect_equal(unique(nodes_df$vertex.color[nodes_df$graphics_type == "line"]), "transparent")
  expect_equal(unique(nodes_df$vertex.color[nodes_df$type == "group"]), "transparent")
})
