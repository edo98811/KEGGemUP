suppressMessages({
  library(KEGGemUP)
  library(igraph)
  library(xml2)
  library(BiocFileCache)
})

# existing nodes: 1111, 2222, 3333, K00001, tst00002, C00001, C00002, C00003
# Existing genes and compounds (subset of the KGML)
nodes_A <- c("1111", "C00001")

# Includes some non-existing IDs (to test missing handling)
nodes_B <- c("3333", "9999", "C00008", "K00001")

# Duplicates within the same vector
nodes_C <- c("1111", "1111", "C00003", "C00002", "C99999")

# All compounds including a non-existent one
nodes_compounds <- c("C00001", "C00002", "C00003", "C00123")


nodes_A_df <- data.frame(
  KEGGID = nodes_A,
  log2FoldChange = rnorm(length(nodes_A), mean = 0, sd = 1)
)
rownames(nodes_A_df) <- nodes_A_df$KEGGID

nodes_B_df <- data.frame(
  KEGG = nodes_B,
  log2FC = rnorm(length(nodes_B), mean = 0, sd = 1)
)
nodes_C_df <- data.frame(
  KEGG_ids = nodes_C,
  log2FoldChange = rnorm(length(nodes_C), mean = 0, sd = 1)
)
nodes_compounds_df <- data.frame(
  KEGG = nodes_compounds,
  log2FC = rnorm(length(nodes_compounds), mean = 0, sd = 1)
)


# Mixed types and random naming
de_results_list <- list(
  transcriptomics = list(
    de_table = nodes_A_df,
    value_column = "log2FoldChange",
    feature_column = "KEGGID"
  ),
  proteomics = list(
    de_table = nodes_C_df,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  ),
  metabolomics = list(
    de_table = nodes_compounds_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  ),
  weird_case = list(
    de_table = nodes_B_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  )
)

# throw_warning <- names(all_de_test_lists)[c(3, 4, 5)]
# expected_warnings <- setNames(c(2, 2, 4), throw_warning)
kgml_path_01 <- system.file("extdata", "test01.xml", package = "KEGGemUP")
kgml_path_real <- system.file("extdata", "hsa04010.xml", package = "KEGGemUP")
kgml_path_02 <- system.file("extdata", "test02.xml", package = "KEGGemUP")
kgml_path_broken <- system.file("extdata", "broken.xml", package = "KEGGemUP")
kgml_path_empty <- system.file("extdata", "empty.xml", package = "KEGGemUP")
kgml_path_no_edges <- system.file("extdata", "no_edges.xml", package = "KEGGemUP")
kgml_path_invalid <- system.file("extdata", "no_kgml.xml", package = "KEGGemUP")

kgml_steps <- readRDS(system.file("extdata", "kgml_parsing_steps.rds", package = "KEGGemUP"))
xml_example <- xml2::read_xml(kgml_path_02)
expected_graphs <- readRDS(system.file("extdata", "kegg_to_graph_expected.rds", package = "KEGGemUP"))
bfc_path <- tools::R_user_dir("BiocFileCache", which = "cache")
bfc <- BiocFileCache(cache = file.path(bfc_path, "test"), ask = FALSE)

vertices_df_basic <- data.frame(
  id = c("n1", "n2", "n3", "n4", "n5", "n6"),
  type = c("gene", "compound", "compound", "compound", "compound", "gene"),
  KEGG = c(NA, "C00001", "C99999", "G00001", "G99999", "00001"),
  label = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
  stringsAsFactors = FALSE
)

# devtools::load_all()

# # Tests for expected outputs of kegg_to_graph and map_results_to_graph
# kgml_steps <- list()

# kgml_steps$nodes <- parse_kgml_nodes(xml_example, kegg_node_defaults())
# kgml_steps$groups <- parse_kgml_groups(xml_example, kegg_node_defaults())
# kgml_steps$line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())

# kgml_steps$all_nodes <- rbind(
#   kgml_steps$nodes,
#   kgml_steps$groups,
#   kgml_steps$line_nodes
# )

# kgml_steps$relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults())
# kgml_steps$reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults())
# kgml_steps$line_edges <- parse_kgml_lines_edges(kgml_steps$line_nodes, kegg_edge_defaults())

# kgml_steps$all_edges <- rbind(
#   kgml_steps$relations_edges,
#   kgml_steps$reactions_edges,
#   kgml_steps$line_edges
# )

# bfc_path <- tools::R_user_dir("BiocFileCache", which = "cache")
# bfc_map <- BiocFileCache(cache = file.path(bfc_path, "mappings"), ask = FALSE)
# g_1 <- build_kegg_graph(kgml_path_01, pathway_name = "hsa00001", bfc_map = bfc)
# g_2 <- build_kegg_graph(kgml_path_02, pathway_name = "hsa00001", bfc_map = bfc)

# kgml_steps$g_test_01 <- g_1
# kgml_steps$g_test_02 <- g_2

# saveRDS(kgml_steps, file = "inst/extdata/kgml_parsing_steps.rds")

# g_test_01 <- kegg_to_graph(
#   pathway_id = "hsa00001",
#   kgml_file = kgml_path_01
# )

# g_test_02 <- kegg_to_graph(
#   pathway_id = "hsa00001",
#   kgml_file = kgml_path_02
# )

# expected <- list(
#   g_test_01 = g_test_01,
#   g_test_02 = g_test_02
# )

# g_test_01_mapped <- map_results_to_graph(
#   g = g_test_01,
#   de_results = de_results_list
# )
# g_test_02_mapped <- map_results_to_graph(
#   g = g_test_02,
#   de_results = de_results_list
# )

# expected$g_test_01_mapped <- g_test_01_mapped
# expected$g_test_02_mapped <- g_test_02_mapped

# vis_graph_01 <- make_kegg_visNetwork(g_test_01_mapped)
# vis_graph_02 <- make_kegg_visNetwork(g_test_02_mapped)

# expected$visNetwork_test_01_mapped <- vis_graph_01
# expected$visNetwork_test_02_mapped <- vis_graph_02

# saveRDS(expected, file = "inst/extdata/kegg_to_graph_expected.rds")



