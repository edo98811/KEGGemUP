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
  mixed_omics           = de_results_list_2, # Mixed omics, column name variations
  duplicates_overlap    = de_results_list_3, # Duplicate / overlapping feature IDs
  all                   = de_results_list_4 # Large mixed test case
)

throw_warning <- names(all_de_test_lists)[c(2, 3, 4)]
expected_warnings <- setNames(c(2, 2, 4), throw_warning)

kgml_path <- system.file("extdata", "test01.xml", package = "KEGGemUP")
kgml_pah_real_example <- system.file("extdata", "hsa04010.xml", package = "KEGGemUP")

# edges_df <- parse_kgml_relations(kgml_path)
# write.table(edges_df, system.file("extdata", "test01.xml_edges.csv", package = "KEGGemUP"),  sep = ";", row.names = FALSE)

# nodes_df <- parse_kgml_entries(kgml_path)
# write.table(nodes_df, system.file("extdata", "test01.xml_nodes.csv", package = "KEGGemUP"),  sep = ";", row.names = FALSE)

nodes_df_path <- system.file("extdata", "test01.xml_nodes.csv", package = "KEGGemUP")
edges_df_path <- system.file("extdata", "test01.xml_edges.csv", package = "KEGGemUP")
kgml_path_empty <- system.file("extdata", "empty_edges.xml", package = "KEGGemUP")

# Expected nodes
expected_nodes <- as.data.frame(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
expected_nodes_cols <- add_columns_nodes_df(expected_nodes)
expected_nodes_cols$KEGG <- vapply(expected_nodes_cols$kegg_name, remove_kegg_prefix_str, FUN.VALUE = character(1))

# Expected edges
expected_edges <- as.data.frame(read.csv(edges_df_path, sep = ";", colClasses = "character"))

# Empty edges
empty_edges <- data.frame(
  from = character(0),
  to = character(0),
  type = character(0),
  subtype = character(0),
  stringsAsFactors = FALSE
)

nodes_df_basic <- data.frame(
  id = c("n1", "n2", "n3", "n4", "n5", "n6"),
  type = c("gene", "compound", "compound", "compound", "compound", "gene"),
  KEGG = c(NA, "C00001", "C99999", "G00001", "G99999", "00001"),
  label = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
  stringsAsFactors = FALSE
)

compounds_file <- system.file(
  "extdata", "compounds.rds",
  package = "KEGGemUP"
)

glycan_file <- system.file(
  "extdata", "glycans.rds",
  package = "KEGGemUP"
)

# Load mapping and mock get_compounds to return it
real_compounds <- readRDS(compounds_file)

real_glycans <- readRDS(glycan_file)