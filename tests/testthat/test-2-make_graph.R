library(mockery)

test_that("expand_keggs handles multiple KEGG IDs correctly", {
  input <- data.frame(
    name = c(1, 2, 3),
    KEGG = c("hsa:1234;hsa:5678", "cpd:C00022", "ko:K00001;ko:K00002;ko:K00003"),
    stringsAsFactors = FALSE
  )
  expected_output <- data.frame(name = c(1, 1, 2, 3, 3, 3), KEGG = c("1234", "5678", "C00022", "K00001", "K00002", "K00003"))
  actual_output <- expand_keggs(input)
  expect_equal(actual_output, expected_output)
})

test_that("remove_kegg_prefix_str removes prefixes and handles multiple IDs", {
  input <- c("hsa:1234 hsa:5678", "cpd:C00022", "ko:K00001 ko:K00002 ko:K00003")
  expected_output <- c("1234;5678", "C00022", "K00001;K00002;K00003")
  actual_output <- vapply(input, remove_kegg_prefix_str, FUN.VALUE = character(1), USE.NAMES = FALSE)
  expect_equal(actual_output, expected_output)
})

test_that("parse_kgml_edges load relationsps correctly", {
  edges_df <- parse_kgml_relations(kgml_path)

  edges_df_expected <- tibble::as_tibble(read.csv(edges_df_path, sep = ";", colClasses = "character"))
  # sapply(colnames(edges_df), function(x) class(edges_df[[x]]))

  expect_equal(edges_df, edges_df_expected)
})

test_that("parse_kgml_entries load nodes correctly", {
  nodes_df <- parse_kgml_entries(kgml_path)

  nodes_df_expected <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))

  expect_equal(nodes_df, nodes_df_expected)
})

# test_that("combine_results_in_dataframe correctly merges DE results", {
#   # Run the function using your already-defined inputs
#   result <- combine_results_in_dataframe(de_results_list)

#   # Check that result is a data frame
#   expect_true(is.data.frame(result))

#   # Check expected columns
#   expect_equal(colnames(result), c("KEGG", "value", "source"))

#   # Check that 'source' column matches names of results_list
#   expect_setequal(unique(result$source), names(de_results_list))

#   # Check KEGG IDs no longer have the prefix (assuming remove_kegg_prefix removes 'path:')
#   expect_false(any(grepl("^hsa:", result$KEGG)))

#   # Check for non-missing values
#   expect_false(any(is.na(result$KEGG)))
#   expect_false(any(is.na(result$value)))

#   # Check that all entries from each source are included
#   expected_nrows <- sum(sapply(de_results_list, function(x) nrow(x$de_table)))
#   expect_equal(nrow(result), expected_nrows)
# })

test_that("combine_results_in_dataframe correctly merges DE results across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    # message(glue::glue("🧩 Testing combine_results_in_dataframe() for scenario: {test_name}"))

    # Run the function for this test list
    result <- combine_results_in_dataframe(de_list)

    # --- Structure checks ---
    expect_true(is.data.frame(result), info = test_name)
    expect_equal(colnames(result), c("KEGG", "value", "source"),
                 info = glue::glue("{test_name}: unexpected column names"))

    # --- Source tracking ---
    expect_setequal(unique(result$source), names(de_list))

    # --- KEGG ID cleanup ---
    expect_false(any(grepl("^hsa:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^cpd:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^path:", result$KEGG)), info = test_name)

    # --- Missing value checks ---
    expect_false(any(is.na(result$KEGG)), info = test_name)
    expect_false(any(is.na(result$value)), info = test_name)

    # --- Row count check ---
    expected_nrows <- sum(vapply(de_list, function(x) nrow(x$de_table), numeric(1)))
    expect_equal(nrow(result), expected_nrows,
                 info = glue::glue("{test_name}: unexpected number of merged rows"))
  })
})


test_that("download_kgml caches and returns file path", {

  # Locate real KGML file in package
  xml_file <- system.file("extdata", "hsa04010.xml", package = "KEGGemUP")

  expect_true(file.exists(xml_file)) # sanity check

  # Read its raw bytes (what get_kgml would normally return)
  xml_content <- readBin(xml_file, what = "raw", n = file.info(xml_file)$size)

  # Create mock get_kgml() returning that real file content
  mock_get <- mock(xml_content)

  # Temporary BiocFileCache location
  bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kgml = mock_get,  # replace only inside this block
    {
      file <- download_kgml("hsa40010", bfc)

      expect_true(file.exists(file))

      # verify mock was called exactly once
      expect_equal(length(mock_args(mock_get)), 1)

      # check that file contains KGML XML
      expect_true(check_valid_kgml(file))
    }
  )
})

test_that("add_compound_names caches and assigns compound names", {

  # Fake nodes df (compound + non-compound)
  nodes_df <- data.frame(
    id = c("n1", "n2", "n3"),
    type = c("gene", "compound", "compound"),
    graphics_name = c(NA, "C00006", "C99999"), 
    label = c("", "", ""),
    stringsAsFactors = FALSE
  )

  # Locate packaged compound mapping
  compounds_file <- system.file(
    "extdata", "compounds.rds",
    package = "KEGGemUP"
  )

  expect_true(file.exists(compounds_file))

  # Load mapping and mock get_compounds to return it
  real_compounds <- readRDS(compounds_file)
  mock_get <- mock(real_compounds)

  # Temporary BiocFileCache
  bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kegg_compounds = mock_get,
    {
      res <- add_compound_names(nodes_df, bfc)
    }
  )
  # ---- assertions ----

  # get_compounds must be called exactly once
  expect_equal(length(mock_args(mock_get)), 1)

  # gene node unchanged
  expect_equal(res$label[1], "")

  # known compound gets correct label
  expect_equal(res$label[2], gsub(";.*", "", as.character(real_compounds[["C00006"]])))

  # unknown compound keeps original ID
  expect_equal(res$label[3], "C99999")
})
