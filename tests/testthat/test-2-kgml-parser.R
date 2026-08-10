test_that("parse_kgml_nodes returns correct nodes", {
  vertices_df <- parse_kgml_nodes(xml_example, kegg_vertex_defaults())
  expect_equal(vertices_df, expected_kgml_steps$nodes)
})

test_that("parse_kgml_groups returns correct group nodes", {
  groups_df <- parse_kgml_groups(xml_example, kegg_vertex_defaults())
  expect_equal(groups_df, expected_kgml_steps$groups)
})

test_that("parse_kgml_lines returns correct line nodes", {
  line_nodes <- parse_kgml_lines(xml_example, kegg_vertex_defaults())
  expect_equal(line_nodes, expected_kgml_steps$line_nodes)
})

test_that("parse_kgml_lines_edges returns correct edges from lines", {
  line_edges <- parse_kgml_lines_edges(expected_kgml_steps$line_nodes, kegg_edge_defaults())
  expect_equal(line_edges, expected_kgml_steps$line_edges)
})

test_that("parse_kgml_relations returns correct edges from relations", {
  relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults())
  expect_equal(relations_edges, expected_kgml_steps$relations_edges)
})

test_that("parse_kgml_reactions returns correct edges from reactions", {
  reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults())
  expect_equal(reactions_edges, expected_kgml_steps$reactions_edges)
})

test_that("complete_kgml_reactions returns correct edges from reactions", {
  all_reaction_edges <- complete_kgml_reactions(expected_kgml_steps$nodes, expected_kgml_steps$reactions_edges, kegg_edge_defaults())
  expect_equal(all_reaction_edges, expected_kgml_steps$completed_reactions)
})

test_that("parse_kgml_reactions returns correct edges from reactions if edges is NULL", {
  reactions_edges <- parse_kgml_reactions(xml_no_edges, kegg_edge_defaults())
  expect_equal(reactions_edges, NULL)
})

test_that("parse_kgml_relations returns correct edges from relations if edges is NULL", {
  relations_edges <- parse_kgml_relations(xml_no_edges, kegg_edge_defaults())
  expect_equal(relations_edges, NULL)
})

test_that("complete_kgml_reactions returns correct edges from reactions if edges is NULL", {
  all_reaction_edges <- complete_kgml_reactions(expected_kgml_steps$nodes, NULL, kegg_edge_defaults())
  expect_equal(all_reaction_edges, NULL)
})

test_that("build_kegg_graph constructs the expected graph (pathway 01)", {
  g <- build_kegg_graph(kgml_path_01, pathway_name = "hsa00001", bfc_map = bfc)
  # expect_true(igraph::identical_graphs(g, expected_kgml_steps$g_test_01))

  expect_equal(
    igraph::graph_attr(g),
    igraph::graph_attr(expected_kgml_steps$g_test_01)
  )

  # Because the vertex order can be different
  v1 <- igraph::as_data_frame(g, what = "vertices")
  v2 <- igraph::as_data_frame(expected_kgml_steps$g_test_01, what = "vertices")

  v1 <- v1[order(v1$name), ]
  v2 <- v2[order(v2$name), ]

  expect_equal(v1, v2)

  e1 <- igraph::as_data_frame(g, what = "edges")
  e2 <- igraph::as_data_frame(expected_kgml_steps$g_test_01, what = "edges")

  e1 <- e1[order(e1$from, e1$to), ]
  e2 <- e2[order(e2$from, e2$to), ]

  expect_equal(e1, e2)
  expect_equal(igraph::graph_attr(g, "title"), "hsa00001")
  # expect_equal(igraph::graph_attr(g, "type"), "KEGG_Pathway")
})

test_that("add_group correctly assigns group labels", {
  nodes <- expected_kgml_steps$all_nodes
  nodes$label <- nodes$name
  nodes_updated <- add_group(nodes)
  group_nodes <- nodes_updated[nodes_updated$type == "group", ]
  expect_false(any(is.na(group_nodes$group)))

  expect_true(all(group_nodes$group == "1, 2"))
  group_idx <- which(nodes_updated$type == "group")
  expect_false(is.na(nodes_updated$x[group_idx]))
  expect_false(is.na(nodes_updated$y[group_idx]))
  expect_true(nodes_updated$x[group_idx] != 0)
  expect_true(nodes_updated$y[group_idx] != 0)
})

test_that("add_labels works correctly", {
  vertices_df <- data.frame(
    KEGG = c("cpd:C00001", "cpd:C00099", "gl:G00001", "ko:K00001", "ko:K99999", "ec:1.1.1.1", "ec:9.9.9.9"),
    ids_for_mapping = c("C00001", "C00099", "G00001", "K00001", "K99999", "1.1.1.1", "9.9.9.9"),
    graphics_name = c("Water", "beta-Alanine, test", "N-Acetyl-D-glucosaminyldiphosphodolichol", "alcohol dehydrogenase", NA, "alcohol dehydrogenase", NA)
  )

  nodes_out <- add_node_labels(vertices_df, bfc = bfc)
  expect_equal(nodes_out$label[1], "H2O")
  expect_equal(nodes_out$label[2], "beta-Alanine")
  expect_equal(nodes_out$label[3], "N-Acetyl-D-glucosaminyldiphosphodolichol")
  # expect_equal(nodes_out$label[4], "E1.1.1.1, adh")
  expect_equal(nodes_out$label[5], "K99999")
  expect_equal(nodes_out$label[6], "alcohol dehydrogenase")
  expect_equal(nodes_out$label[7], "9.9.9.9")
})

test_that("add_reaction_labels works correctly", {
  vertices_df <- data.frame(
    reaction = c("rn:R00001", "rn:R00099", NA, "rn:R00002")
  )

  nodes_out <- add_reaction_labels(vertices_df, bfc = bfc)
  expect_equal(nodes_out$reaction_label[1], "polyphosphate polyphosphohydrolase")
  expect_equal(nodes_out$reaction_label[2], "Cob(I)alamin <=> Cob(II)alamin")
  expect_true(is.na(nodes_out$reaction_label[3]))
  expect_equal(nodes_out$reaction_label[4], "reduced ferredoxin:dinitrogen oxidoreductase (ATP-hydrolysing)")

  expect_equal(nodes_out$reaction_link[1], "https://www.kegg.jp/dbget-bin/www_bget?R00001")
  expect_equal(nodes_out$reaction_link[2], "https://www.kegg.jp/dbget-bin/www_bget?R00099")
  expect_true(is.na(nodes_out$reaction_link[3]))
  expect_equal(nodes_out$reaction_link[4], "https://www.kegg.jp/dbget-bin/www_bget?R00002")
})
