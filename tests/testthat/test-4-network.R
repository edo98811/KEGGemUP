
# # The functions to buil the graph can handle empty or malformed KGML files, 
# # and will provide informative error messages or warnings when such cases are encountered.
# test_that("build_kegg_graph handles empty KGML files", {
#   expect_error(
#     build_kegg_graph(kgml_path_broken, "Broken_Pathway", bfc = bfc),
#     regexp = "Failed to read XML file"
#   )
# })
# test_that("build_kegg_graph handles empty KGML files", {
#   expect_error(
#     build_kegg_graph(kgml_path_empty, "Empty_Pathway", bfc = bfc),
#     regexp = "Failed to read XML file"
#   )
# })
# test_that("build_kegg_graph handles empty KGML files", {
#   expect_warning(
#     build_kegg_graph(kgml_path_no_edges, "No_edges_Pathway", bfc = bfc),
#     regexp = "No edges in graph"
#   )
# })


# test_that("build_kegg_graph works correctly", {
#   g <- build_kegg_graph(kgml_path_01, "hsa00001", bfc = bfc)

#   expect_true(inherits(g, "igraph"))
#   vertices_df <- as_data_frame(g, what = "vertices")
#   edges_df <- as_data_frame(g, what = "edges")

#   expect_true(all(edges_df$from %in% vertices_df$name))
#   expect_true(all(edges_df$to %in% vertices_df$name))

#   expect_equal(igraph::graph_attr(g, "title"), "hsa00001")
#   expect_equal(igraph::graph_attr(g, "type"), "KEGG_Pathway")

#   expect_true(all(is.character(vertices_df$ids_for_mapping)))
#   expect_true(all(vertices_df$ids_for_mapping == "" | nchar(vertices_df$ids_for_mapping)))
# })

# test_that("make_igraph_graph handles edge cases", {

#   ## Empty edges (to be sure)
#   vertices_df <- data.frame(name = c("1", "2"), label = c("A", "B"))
#   edges_df <- data.frame(from = character(), to = character())

#   expect_warning(
#     g <- make_igraph_graph(vertices_df, edges_df, "hsa00001"),
#     regexp = "No edges in graph"
#   )

#   expect_true(inherits(g, "igraph"))
#   expect_equal(igraph::vcount(g), 2)
#   expect_equal(igraph::ecount(g), 0)

#   # Same with NULL edges (what is expected)
#   expect_warning(
#     g <- make_igraph_graph(vertices_df, NULL, "hsa00001"),
#     regexp = "No edges in graph"
#   )

#   expect_true(inherits(g, "igraph"))
#   expect_equal(igraph::vcount(g), 2)
#   expect_equal(igraph::ecount(g), 0)

#   vertices_df <- data.frame(name = c("3", "1", "2"), label = c("Z", "A", "B"))
#   edges_df <- data.frame(from = c("3", "1"), to = c("1", "2"))

#   g <- make_igraph_graph(vertices_df, edges_df, "hsa00001")

#   vertex_labels <- igraph::V(g)$label
#   expect_equal(vertex_labels, sort(vertex_labels))
# })

# test_that("make_igraph_graph works correctly", {

#   # Normal case
#   vertices_df <- data.frame(name = c("1", "2", "3"), label = c("A", "B", "C"))
#   edges_df <- data.frame(from = c("1", "2"), to = c("2", "3"))

#   g <- make_igraph_graph(vertices_df, edges_df, "hsa00001")

#   expect_true(inherits(g, "igraph"))
#   expect_equal(igraph::vcount(g), 3)
#   expect_equal(igraph::ecount(g), 2)
#   })


  