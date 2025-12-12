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
  edges_df <- suppressMessages(parse_kgml_relations(kgml_path))

  expect_equal(edges_df, expected_edges)
})

test_that("parse_kgml_entries loads empty edges  correctly", {
  expect_warning(edges_df <- suppressMessages(parse_kgml_relations(kgml_path_empty)), "No relations found in KGML file.")

  expect_true(nrow(edges_df) == 0)
  expect_true(inherits(edges_df, "data.frame"))
  expect_equal(colnames(edges_df), c("from", "to", "type", "subtype", "rel_value"))
})

test_that("parse_kgml_entries load nodes correctly", {
  nodes_df <- suppressMessages(parse_kgml_entries(kgml_path))

  expect_equal(nodes_df, expected_nodes)
})

test_that("combine_results_in_dataframe correctly merges DE results across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    # Run the function for this test list
    result <- combine_results_in_dataframe(de_list)

    # Structure checks
    expect_true(is.data.frame(result), info = test_name)
    expect_equal(colnames(result), c("KEGG", "plot_value", "source"),
      info = paste0(test_name, "unexpected column names")
    )

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
      info = paste0(test_name, "unexpected number of merged rows")
    )
  })
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


# Example test using mockery
test_that("download_kgml works without hitting KEGG API", {
  
  # Create temporary BiocFileCache and directory
  bfc <- BiocFileCache(tempdir(), ask = FALSE)
  temp_dir <- tempdir()
  valid_pathway <- "hsa04110"

  # Create a fake XML response
  fake_xml <- read_xml("<pathway name='fake_pathway'></pathway>")

  # Create mocks for httr2 functions
  mock_req_retry <- function(req, max_tries = 1) req
  mock_req_perform <- function(req, error_call = FALSE) structure(list(), class = "response")
  mock_resp_is_error <- function(resp) FALSE
  mock_resp_body_xml <- function(resp) fake_xml
  mock_resp_status <- function(resp) 200

  # Use with_mocked_bindings to replace httr2 functions within this scope
  kgml_path <- suppressMessages(with_mocked_bindings(
    download_kgml(
      valid_pathway,
      directory = temp_dir
    ),
    request = function(url) list(url = url),
    req_retry = mock_req_retry,
    req_perform = mock_req_perform,
    resp_is_error = mock_resp_is_error,
    resp_body_xml = mock_resp_body_xml,
    resp_status = mock_resp_status
  ))

  # Directory mode checks
  expect_true(file.exists(kgml_path))
  expect_true(grepl("\\.xml$", kgml_path))

  # Cache mode with same mocks
  kgml_cache <- suppressMessages(with_mocked_bindings(
    download_kgml(
      valid_pathway,
      bfc = bfc
    ),
    request = function(url) list(url = url),
    req_retry = mock_req_retry,
    req_perform = mock_req_perform,
    resp_is_error = mock_resp_is_error,
    resp_body_xml = mock_resp_body_xml,
    resp_status = mock_resp_status
  ))

  expect_true(file.exists(kgml_cache))
  expect_true(grepl("\\.xml$", kgml_cache))

  # Second download (should reuse cache)
  expect_message(
    kgml_cache2 <- with_mocked_bindings(
      download_kgml(valid_pathway, bfc = bfc),
      request = function(url) list(url = url),
      req_retry = mock_req_retry,
      req_perform = mock_req_perform,
      resp_is_error = mock_resp_is_error,
      resp_body_xml = mock_resp_body_xml,
      resp_status = mock_resp_status
    ),
    "Using cached KEGG KGML"
  )

  # Test that cached path is reused
  expect_equal(kgml_cache, kgml_cache2)
})
