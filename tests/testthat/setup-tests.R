# existing nodes: 1111, 2222, 3333, K00001, tst00002, C00001, C00002, C00003
# Existing genes and compounds (subset of the KGML)
nodes_A <- c("1111", "C00001")

# Includes some non-existing IDs (to test missing handling)
nodes_B <- c("3333", "9999", "C00008", "K00001")

# Duplicates within the same vector
nodes_C <- c("1111", "1111", "C00003", "C00002", "C99999")

# Mix of existing and new, with overlap across vectors
nodes_D <- c("C00003", "C00008", "3333", "tst00002")

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
nodes_D_df <- data.frame(
  KEGG_ids = nodes_D,
  log2FoldChange = rnorm(length(nodes_D), mean = 0, sd = 1)
)
nodes_compounds_df <- data.frame(
  KEGG = nodes_compounds,
  log2FC = rnorm(length(nodes_compounds), mean = 0, sd = 1)
)

# --- EXAMPLE FAKE DE RESULTS LISTS ---

# Basic case — two datasets, consistent and simple
de_results_list_1 <- list(
  genes = list(
    de_table = nodes_A_df,
    value_column = "log2FoldChange",
    feature_column = "KEGGID"
  ),
  metabolites = list(
    de_table = nodes_B_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  )
)

de_results_list_rownames <- list(
  genes = list(
    de_table = nodes_A_df,
    value_column = "log2FoldChange",
    feature_column = "rownames"
  ),
  metabolites = list(
    de_table = nodes_B_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  )
)


# Mixed column names and redundant identifiers
de_results_list_2 <- list(
  transcr = list(
    de_table = nodes_C_df,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  ),
  proteins = list(
    de_table = nodes_B_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  ),
  metabolome = list(
    de_table = nodes_compounds_df,
    value_column = "log2FC",
    feature_column = "KEGG"
  )
)

# Duplicates and cross-referenced names 
de_results_list_3 <- list(
  group1 = list(
    de_table = nodes_D_df,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  ),
  group2 = list(
    de_table = nodes_C_df,
    value_column = "log2FoldChange",
    feature_column = "KEGG_ids"
  )
)

# Mixed types and random naming 
de_results_list_4 <- list(
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

all_de_test_lists <- list(
  genes_metabolites     = de_results_list_1, # Basic and consistent
  using_rownames        = de_results_list_rownames, # Basic and consistent with rownames as feature_column
  mixed_omics           = de_results_list_2, # Mixed omics, column name variations
  duplicates_overlap    = de_results_list_3, # Duplicate / overlapping feature IDs
  all                   = de_results_list_4 # Large mixed test case
)

throw_warning <- names(all_de_test_lists)[c(3, 4, 5)]
expected_warnings <- setNames(c(2, 2, 4), throw_warning)

kgml_path <- system.file("extdata", "test01.xml", package = "KEGGemUP")
ref_graph <- readRDS(system.file("extdata", "test01_reference_graph.rds", package = "KEGGemUP"))
xml_example <- xml2::read_xml(kgml_path)
kgml_processing_steps <- readRDS(system.file("extdata", "kgml_parsing_steps.rds", package = "KEGGemUP"))

# Empty edges
empty_edges <- data.frame(
  from = character(0),
  to = character(0),
  type = character(0),
  relation_subtype = character(0),
  stringsAsFactors = FALSE
)

nodes_df_basic <- data.frame(
  id = c("n1", "n2", "n3", "n4", "n5", "n6"),
  type = c("gene", "compound", "compound", "compound", "compound", "gene"),
  KEGG = c(NA, "C00001", "C99999", "G00001", "G99999", "00001"),
  label = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
  stringsAsFactors = FALSE
)

# make working nodes
kgml_steps <- list()

kgml_steps$nodes <- parse_kgml_nodes(xml_example, kegg_node_defaults())
kgml_steps$groups <- parse_kgml_groups(xml_example, kegg_node_defaults())
kgml_steps$line_nodes <- parse_kgml_lines(xml_example, kegg_node_defaults())

kgml_steps$all_nodes <- rbind(
  kgml_steps$nodes,
  kgml_steps$groups,
  kgml_steps$line_nodes
)

kgml_steps$relations_edges <- parse_kgml_relations(xml_example, kegg_edge_defaults())
kgml_steps$reactions_edges <- parse_kgml_reactions(xml_example, kegg_edge_defaults())
kgml_steps$line_edges <- parse_kgml_lines_edges(kgml_steps$line_nodes, kegg_edge_defaults())

kgml_steps$all_edges <- rbind(
  kgml_steps$relations_edges,
  kgml_steps$reactions_edges,
  kgml_steps$line_edges
)

saveRDS(kgml_steps, file = "kgml_parsing_steps.rds")


kgml_path <- system.file("extdata", "test01.xml", package = "YourPackageName")

g_test <- kegg_to_graph(
  pathway_id = "hsa:TEST01",
  scaling_factor = 1.0,
  simplified_graph = TRUE,
  kgml_file = kgml_path
)

v_test <- plot_visNetwork_kegg(g_mapped)
