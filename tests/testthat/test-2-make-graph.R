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

test_that("parse_kgml_nodes + groups + lines load nodes correctly", {
  expected_nodes <- readRDS(system.file("extdata", "test01.xml_nodes_combined.rds", package = "KEGGemUP"))

  # Parse nodes
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml_example, kegg_node_defaults()))
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())
  nodes_df <- rbind(nodes_df, line_nodes)

  expect_equal(nodes_df, expected_nodes)
})

test_that("parse_kgml_relations + reactions + line edges load edges correctly", {
  expected_edges <- readRDS(system.file("extdata", "test01.xml_example_edges_combined.rds", package = "KEGGemUP"))

  # Parse edges
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())
  edges_df <- parse_kgml_relations(xml_example, kegg_edge_defaults())
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml_example, kegg_edge_defaults()))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults()))

  expect_equal(edges_df, expected_edges)
})

test_that("build_kegg_graph constructs the expected graph", {
  bfc <- BiocFileCache(tempfile(), ask = FALSE)
  g <- build_kegg_graph(kgml_path, pathway_name = "Test_Pathway", bfc_map = bfc)

  # Compare igraph objects
  # We use igraph::identical_graphs to ensure topology + attributes match
  expect_true(igraph::identical_graphs(g, expected_graph))

  # Optionally, check key graph attributes
  expect_equal(igraph::graph_attr(g, "title"), "Test_Pathway")
  expect_equal(igraph::graph_attr(g, "type"), "Test_Pathway")
})
