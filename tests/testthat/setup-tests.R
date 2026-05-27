# loading required packages for setup -------------------------------------

suppressMessages({
  library(KEGGemUP)
  library(igraph)
  library(xml2)
  library(BiocFileCache)
})

# definition of expected files and co -------------------------------------

kgml_path_01 <- system.file("extdata", "test01.xml", package = "KEGGemUP")
kgml_path_broken <- system.file("extdata", "broken.xml", package = "KEGGemUP")
kgml_path_empty <- system.file("extdata", "empty.xml", package = "KEGGemUP")
kgml_path_no_edges <- system.file("extdata", "no_edges.xml", package = "KEGGemUP")

xml_example <- xml2::read_xml(kgml_path_01)
xml_no_edges <- xml2::read_xml(kgml_path_no_edges)

bfc_path <- tools::R_user_dir("BiocFileCache", which = "cache")
bfc <- BiocFileCache(cache = file.path(bfc_path, "test"), ask = FALSE)




# specifying some nodes ---------------------------------------------------

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


# DE-like content ---------------------------------------------------------

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



# creating the rds objects for the testing --------------------------------

## Tests for expected outputs of create_kegg_graph and map_results_to_graph

expected_kgml_steps <- list()

expected_kgml_steps$nodes <- parse_kgml_nodes(xml_example, kegg_vertex_defaults())
expected_kgml_steps$groups <- parse_kgml_groups(xml_example, kegg_vertex_defaults())
expected_kgml_steps$line_nodes <- parse_kgml_lines(xml_example, kegg_vertex_defaults())

expected_kgml_steps$all_nodes <- rbind(
  expected_kgml_steps$nodes,
  expected_kgml_steps$groups,
  expected_kgml_steps$line_nodes
)

expected_kgml_steps$relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults())
expected_kgml_steps$reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults())
expected_kgml_steps$line_edges <- parse_kgml_lines_edges(expected_kgml_steps$line_nodes,
                                                kegg_edge_defaults())
expected_kgml_steps$completed_reactions <-
  complete_kgml_reactions(expected_kgml_steps$nodes, expected_kgml_steps$reactions_edges, kegg_edge_defaults())

expected_kgml_steps$all_edges <- rbind(
  expected_kgml_steps$relations_edges,
  expected_kgml_steps$reactions_edges,
  expected_kgml_steps$line_edges,
  complete_kgml_reactions(expected_kgml_steps$nodes, expected_kgml_steps$reactions_edges, kegg_edge_defaults())
)

bfc_path <- tools::R_user_dir("BiocFileCache", which = "cache")
bfc_map <- BiocFileCache(cache = file.path(bfc_path, "mappings"), ask = FALSE)

g_1 <- suppressWarnings({
  build_kegg_graph(kgml_path_01, pathway_name = "hsa00001", bfc_map = bfc_map)
})

expected_kgml_steps$g_test_01 <- g_1

# saveRDS(expected_kgml_steps, file = "inst/extdata/expected_kgml_steps.RDS")

g_test_01 <- suppressWarnings({
  create_kegg_graph(
    pathway_id = "hsa00001",
    kgml_file = kgml_path_01
  )
})

expected_output <- list(
  g_test_01 = g_test_01
)

g_test_01_mapped <- suppressWarnings({
  map_results_to_graph(
    g = g_test_01,
    de_results = de_results_list
  )
})

expected_output$g_test_01_mapped <- g_test_01_mapped

vis_graph_01 <- render_kegg_graph(g_test_01_mapped)

expected_output$visNetwork_test_01_mapped <- vis_graph_01

# saveRDS(expected_output, file = "inst/extdata/expected_output.RDS")

# loading expected parsing steps and outputs - see above for the creation --------

# expected_kgml_steps <- readRDS(
#   system.file("extdata", "expected_kgml_steps.RDS", package = "KEGGemUP")
# )
# expected_output <- readRDS(
#   system.file("extdata", "expected_output.RDS", package = "KEGGemUP")
# )


