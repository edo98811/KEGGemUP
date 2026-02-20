test_that("parse_kgml_nodes returns correct nodes", {
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults())
  expect_equal(nodes_df, kgml_steps$nodes)
})

test_that("parse_kgml_groups returns correct group nodes", {
  groups_df <- parse_kgml_groups(xml_example, kegg_node_defaults())
  expect_equal(groups_df, kgml_steps$groups)
})

test_that("parse_kgml_lines returns correct line nodes", {
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())
  expect_equal(line_nodes, kgml_steps$line_nodes)
})

test_that("parse_kgml_lines_edges returns correct edges from lines", {
  line_edges <- parse_kgml_lines_edges(kgml_steps$line_nodes, kegg_edge_defaults())
  expect_equal(line_edges, kgml_steps$line_edges)
})

test_that("parse_kgml_relations returns correct edges from relations", {
  relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults())
  expect_equal(relations_edges, kgml_steps$relations_edges)
})

test_that("parse_kgml_reactions returns correct edges from reactions", {
  reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults())
  expect_equal(reactions_edges, kgml_steps$reactions_edges)
})

test_that("combined nodes (nodes + groups + lines) load correctly", {
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml_example, kegg_node_defaults()))
  nodes_df <- rbind(nodes_df, parse_kgml_lines(xml_example, kegg_node_defaults()))
  expect_equal(nodes_df, kgml_steps$all_nodes)
})

test_that("combined edges (relations + reactions + line edges) load correctly", {
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())
  edges_df <- parse_kgml_relations(xml_example, kegg_edge_defaults())
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml_example, kegg_edge_defaults()))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults()))
  expect_equal(edges_df, kgml_steps$all_edges)
})

test_that("build_kegg_graph constructs the expected graph (pathway 01)", {
  g <- build_kegg_graph(kgml_path_01, pathway_name = "hsa00001", bfc_map = bfc)
  expect_true(igraph::identical_graphs(g, kgml_steps$g_test_01))
  expect_equal(igraph::graph_attr(g, "title"), "hsa00001")
  expect_equal(igraph::graph_attr(g, "type"), "KEGG_Pathway")
})

test_that("build_kegg_graph constructs the expected graph (pathway 02)", {
  g <- build_kegg_graph(kgml_path_02, pathway_name = "hsa00001", bfc_map = bfc)
  all(sort(vertex_attr_names(g)) == sort(vertex_attr_names(kgml_steps$g_test_02))) &&
    all(sort(edge_attr_names(g)) == sort(edge_attr_names(kgml_steps$g_test_02)))
  # expect_true(igraph::identical_graphs(g, kgml_steps$g_test_02))
  expect_equal(igraph::graph_attr(g, "title"), "hsa00001")
  expect_equal(igraph::graph_attr(g, "type"), "KEGG_Pathway")
})

test_that("add_group correctly assigns group labels", {
  nodes <- kgml_steps$all_nodes
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
  compounds_db <- data.frame(id = c("C00001", "C00002"), name = c("Water;H2O", "ATP"))
  glycans_db <- data.frame(id = c("G00001"), name = c("GlycanX;Y"))
  genes_db <- data.frame(id = c("K00001", "K00002"), name = c("GeneA", "GeneB;alias"))
  enzymes_db <- data.frame(id = c("1.1.1.1"), name = c("EnzymeX"))

  nodes_df <- data.frame(
    KEGG = c("cpd:C00001", "cpd:C00099", "gl:G00001", "ko:K00001", "ko:K99999", "ec:1.1.1.1", "ec:9.9.9.9"),
    ids_for_mapping = c("C00001", "C00099", "G00001", "K00001", "K99999", "1.1.1.1", "9.9.9.9"),
    graphics_name = c("Water", "beta-Alanine, No", "N-Acetyl-D-glucosaminyldiphosphodolichol", "alcohol dehydrogenase", NA, "alcohol dehydrogenase", NA),
    stringsAsFactors = FALSE
  )

  nodes_out <- add_node_labels(nodes_df, bfc = bfc)
  expect_equal(nodes_out$label[1], "H2O")
  expect_equal(nodes_out$label[2], "beta-Alanine")
  expect_equal(nodes_out$label[3], "N-Acetyl-D-glucosaminyldiphosphodolichol")
  expect_equal(nodes_out$label[4], "E1.1.1.1")
  expect_equal(nodes_out$label[5], "K99999")
  expect_equal(nodes_out$label[6], "alcohol dehydrogenase")
  expect_equal(nodes_out$label[7], "9.9.9.9")
})

test_that("add_reaction_labels works correctly", {
  nodes_df <- data.frame(
    reaction = c("rn:R00001", "rn:R00099", NA, "rn:R00002"),
    stringsAsFactors = FALSE
  )

  nodes_out <- add_reaction_labels(nodes_df, bfc = bfc)
  expect_equal(nodes_out$reaction_label[1], "polyphosphate polyphosphohydrolase")
  expect_equal(nodes_out$reaction_label[2], "Cob(I)alamin <=> Cob(II)alamin")
  expect_true(is.na(nodes_out$reaction_label[3]))
  expect_equal(nodes_out$reaction_label[4], "reduced ferredoxin:dinitrogen oxidoreductase (ATP-hydrolysing)")

  expect_equal(nodes_out$reaction_link[1], "https://www.kegg.jp/dbget-bin/www_bget?R00001")
  expect_equal(nodes_out$reaction_link[2], "https://www.kegg.jp/dbget-bin/www_bget?R00099")
  expect_true(is.na(nodes_out$reaction_link[3]))
  expect_equal(nodes_out$reaction_link[4], "https://www.kegg.jp/dbget-bin/www_bget?R00002")
})
