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

# ✅ Basic case — two datasets, consistent and simple
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

# ⚠️ Mixed column names and redundant identifiers
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

# 🔁 Duplicates and cross-referenced names (e.g., same nodes across lists)
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

# 🧪 Mixed types and random naming — stress test
de_results_list_5 <- list(
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
  full_stress           = de_results_list_5  # Large mixed test case 
)

throw_warning <- names(all_de_test_lists)[c(2,3,4)]
expected_warnings <- setNames(c(2, 2, 4), throw_warning)

kgml_path <- system.file("extdata", "test01.xml", package = "KEGGemUP")

edges_df <- parse_kgml_relations(kgml_path)
write.table(edges_df, system.file("extdata", "test01.xml_edges.csv", package = "KEGGemUP"),  sep = ";", row.names = FALSE)

nodes_df <- parse_kgml_entries(kgml_path)
write.table(nodes_df, system.file("extdata", "test01.xml_nodes.csv", package = "KEGGemUP"),  sep = ";", row.names = FALSE)

nodes_df_path <- system.file("extdata", "test01.xml_nodes.csv", package = "KEGGemUP")
edges_df_path <- system.file("extdata", "test01.xml_edges.csv", package = "KEGGemUP")
