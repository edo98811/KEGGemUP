test_that("parse_kgml_nodes returns correct nodes", {
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults)
  expect_equal(nodes_df, kgml_steps$nodes)
})

test_that("parse_kgml_groups returns correct group nodes", {
  groups_df <- parse_kgml_groups(xml_example, kegg_node_defaults)
  expect_equal(groups_df, kgml_steps$groups)
})

test_that("parse_kgml_lines returns correct line nodes", {
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults)
  expect_equal(line_nodes, kgml_steps$line_nodes)
})

test_that("parse_kgml_lines_edges returns correct edges from lines", {
  line_edges <- parse_kgml_lines_edges(kgml_steps$line_nodes, kegg_edge_defaults)
  expect_equal(line_edges, kgml_steps$line_edges)
})

test_that("parse_kgml_relations returns correct edges from relations", {
  relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults)
  expect_equal(relations_edges, kgml_steps$relations_edges)
})

test_that("parse_kgml_reactions returns correct edges from reactions", {
  reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults)
  expect_equal(reactions_edges, kgml_steps$reactions_edges)
})

test_that("combined nodes (nodes + groups + lines) load correctly", {
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults)
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml_example, kegg_node_defaults))
  nodes_df <- rbind(nodes_df, parse_kgml_lines(xml_example, kegg_node_defaults))
  expect_equal(nodes_df, kgml_steps$all_nodes)
})

test_that("combined edges (relations + reactions + line edges) load correctly", {
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults)
  edges_df <- parse_kgml_relations(xml_example, kegg_edge_defaults)
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml_example, kegg_edge_defaults))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults))
  expect_equal(edges_df, kgml_steps$all_edges)
})

test_that("parse_kgml_nodes + groups + lines load nodes correctly", {
  expected_nodes <- readRDS(system.file("extdata", "test01.xml_nodes_combined.rds", package = "KEGGemUP"))

  # Parse nodes
  nodes_df <- parse_kgml_nodes(xml_example, kegg_node_defaults)
  nodes_df <- rbind(nodes_df, parse_kgml_groups(xml_example, kegg_node_defaults))
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults)
  nodes_df <- rbind(nodes_df, line_nodes)

  expect_equal(nodes_df, expected_nodes)
})

test_that("parse_kgml_relations + reactions + line edges load edges correctly", {
  expected_edges <- readRDS(system.file("extdata", "test01.xml_example_edges_combined.rds", package = "KEGGemUP"))

  # Parse edges
  line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults)
  edges_df <- parse_kgml_relations(xml_example, kegg_edge_defaults)
  edges_df <- rbind(edges_df, parse_kgml_reactions(xml_example, kegg_edge_defaults))
  edges_df <- rbind(edges_df, parse_kgml_lines_edges(line_nodes, kegg_edge_defaults))

  expect_equal(edges_df, expected_edges)
})

test_that("build_kegg_graph constructs the expected graph", {
  g <- build_kegg_graph(kgml_path_01, pathway_name = "Test_Pathway", bfc_map = bfc)

  # Compare igraph objects
  # We use igraph::identical_graphs to ensure topology + attributes match
  expect_true(igraph::identical_graphs(g, expected_graph))

  # Optionally, check key graph attributes
  expect_equal(igraph::graph_attr(g, "title"), "Test_Pathway")
  expect_equal(igraph::graph_attr(g, "type"), "Test_Pathway")
})


test_that("add_group correctly assigns group labels", {
  nodes <- kgml_steps$all_nodes

  nodes$label <- nodes$name

  # Run add_group
  nodes_updated <- add_group(nodes)

  # Identify group nodes
  group_nodes <- nodes_updated[nodes_updated$type == "group", ]

  # All group nodes should have non-NA group label
  expect_false(any(is.na(group_nodes$group)))

  # Components of each group node should have the same group label
  for (i in seq_len(nrow(group_nodes))) {
    group_node <- group_nodes[i, ]
    comps <- strsplit(group_node$components, ";", fixed = TRUE)[[1]]
    node_idx <- match(c(comps, group_node$name), nodes_updated$name)
    expect_true(all(nodes_updated$group[node_idx] == group_node$group))
  }

  # Check that non-NA groups have the expected value
  expect_true(all(group_nodes$group == "1;2"))

  group_idx <- which(nodes_updated$type == "group")
  expect_false(is.na(nodes_updated$x[group_idx]))
  expect_false(is.na(nodes_updated$y[group_idx]))
  expect_true(nodes_updated$x[group_idx] != 0)
  expect_true(nodes_updated$y[group_idx] != 0)
})

test_that("build_kegg_graph constructs the expected graph", {
  g <- build_kegg_graph(kgml_path_02, pathway_name = "Test_Pathway", bfc_map = bfc)

  # Compare igraph objects
  # We use igraph::identical_graphs to ensure topology + attributes match
  expect_true(igraph::identical_graphs(g, expected_graph))

  # Optionally, check key graph attributes
  expect_equal(igraph::graph_attr(g, "title"), "Test_Pathway")
  expect_equal(igraph::graph_attr(g, "type"), "Test_Pathway")
})

test_that("add_labels works correctly", {
  compounds_db <- data.frame(id = c("C00001", "C00002"), name = c("Water;H2O", "ATP"))
  glycans_db <- data.frame(id = c("G00001"), name = c("GlycanX;Y"))
  genes_db <- data.frame(id = c("K00001", "K00002"), name = c("GeneA", "GeneB;alias"))
  enzymes_db <- data.frame(id = c("1.1.1.1"), name = c("EnzymeX"))

  nodes_df <- data.frame(
    KEGG = c("cpd:C00001", "cpd:C00099", "cpd:G00001", "ko:K00001", "ko:K99999", "ec:1.1.1.1", "ec:9.9.9.9"),
    ids_for_mapping = c("C00001", "C00099", "G00001", "K00001", "K99999", "1.1.1.1", "9.9.9.9"),
    stringsAsFactors = FALSE
  )

  nodes_out <- add_node_labels(nodes_df, bfc = bfc)

  # Assertions
  expect_equal(nodes_out$label[1], "H2O") # compound found, truncated
  expect_equal(nodes_out$label[2], "beta-Alanine") # compound not found, fallback
  expect_equal(nodes_out$label[3], "N-Acetyl-D-glucosaminyldiphosphodolichol") # glycan, full name
  expect_equal(nodes_out$label[4], "E1.1.1.1, adh") # gene, found
  expect_equal(nodes_out$label[5], "K99999") # gene, not found
  expect_equal(nodes_out$label[6], "alcohol dehydrogenase") # enzyme, found
  expect_equal(nodes_out$label[7], "9.9.9.9") # enzyme, not found
})

test_that("add_reaction_labels works correctly", {

  nodes_df <- data.frame(
    reaction = c("rn:R00001", "rn:R00099", NA, "rn:R00002"),
    stringsAsFactors = FALSE
  )
  devtools::load_all()
  nodes_out <- add_reaction_labels(nodes_df, bfc = bfc)

  # Assertions
  expect_equal(nodes_out$reaction_label[1], "polyphosphate polyphosphohydrolase") # found in DB
  expect_equal(nodes_out$reaction_label[2], "Cob(I)alamin <=> Cob(II)alamin") # not found → fallback
  expect_true(is.na(nodes_out$reaction_label[3])) # NA reaction
  expect_equal(nodes_out$reaction_label[4], "reduced ferredoxin:dinitrogen oxidoreductase (ATP-hydrolysing)") # found

  # Reaction links
  expect_equal(nodes_out$reaction_link[1], "https://www.kegg.jp/dbget-bin/www_bget?R00001")
  expect_equal(nodes_out$reaction_link[2], "https://www.kegg.jp/dbget-bin/www_bget?R00099")
  expect_true(is.na(nodes_out$reaction_link[3]))
  expect_equal(nodes_out$reaction_link[4], "https://www.kegg.jp/dbget-bin/www_bget?R00002")
})
