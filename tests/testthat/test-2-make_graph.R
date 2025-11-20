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

  expect_equal(edges_df, expected_edges)
})

test_that("parse_kgml_entries load nodes correctly", {
  nodes_df <- parse_kgml_entries(kgml_path)

  expect_equal(nodes_df, expected_nodes)
})

test_that("combine_results_in_dataframe correctly merges DE results across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {

    # Run the function for this test list
    result <- combine_results_in_dataframe(de_list)

    # Structure checks
    expect_true(is.data.frame(result), info = test_name)
    expect_equal(colnames(result), c("KEGG", "value", "source"),
                 info = paste0(test_name, "unexpected column names"))

    # Source tracking
    expect_setequal(unique(result$source), names(de_list))

    # KEGG ID cleanup
    expect_false(any(grepl("^hsa:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^cpd:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^path:", result$KEGG)), info = test_name)

    # Missing value checks
    expect_false(any(is.na(result$KEGG)), info = test_name)
    expect_false(any(is.na(result$plot_value)), info = test_name)

    # Row count check
    expected_nrows <- sum(vapply(de_list, function(x) nrow(x$de_table), numeric(1)))
    expect_equal(nrow(result), expected_nrows,
                 info = paste0(test_name, "unexpected number of merged rows"))
  })
})


test_that("download_kgml caches and returns file path", {

  expect_true(file.exists(kgml_pah_real_example))

  # Read its raw bytes (what get_kgml would normally return)
  xml_content <- readBin(kgml_pah_real_example, what = "raw", n = file.info(kgml_pah_real_example)$size)

  mock_get <- mock(xml_content)

  bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kgml = mock_get,  
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

test_that("add_compound_names caches and assigns glycan and compounds names", {

  bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kegg_compounds = mock(real_compounds),
    get_kegg_glycans = mock(real_glycans),
    {
      res <- add_compound_names(nodes_df_basic, bfc)
    }
  )
  
  # gene node unchanged
  expect_equal(res$label[1], NA_character_)

  # known compound gets correct label
  expect_equal(res$label[2], gsub(";.*", "", as.character(real_compounds[["C00001"]])))
  expect_equal(res$label[4], gsub(";.*", "", as.character(real_glycans[["G00001"]])))

  # unknown compound keeps original ID
  expect_equal(res$label[3], "C99999")
})
