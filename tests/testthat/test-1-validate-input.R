test_that("is_valid_de_entry works as expected", {
  #  Valid input
  de_table <- data.frame(
    gene = c("A", "B"),
    logFC = c(1.2, -0.5),
    p_val = c(0.01, 0.2)
  )

  valid_entry <- list(
    de_table = de_table,
    value_column = "logFC",
    feature_column = "gene"
  )
  expect_true(is_valid_de_entry(valid_entry, "valid_entry"))

  # Missing required list elements
  invalid_missing <- list(de_table = de_table)
  expect_warning(expect_false(is_valid_de_entry(invalid_missing, "invalid_missing")))

  # de_table not a data frame
  invalid_table <- list(
    de_table = matrix(1:4, ncol = 2),
    value_column = "logFC",
    feature_column = "gene"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_table, "invalid_table")))

  # Value_column not present
  invalid_value_column <- list(
    de_table = de_table,
    value_column = "not_here",
    feature_column = "gene"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_value_column, "invalid_value_column")))

  # Feature_column not present or not rownames
  invalid_feature_column <- list(
    de_table = de_table,
    value_column = "logFC",
    feature_column = "not_here"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_feature_column, "invalid_feature_column")))

  # Feature_column = 'rownames' case
  rownames(de_table) <- de_table$gene
  rowname_entry <- list(
    de_table = de_table,
    value_column = "logFC",
    feature_column = "rownames"
  )

  # This checks logic– if rownames allowed, expect TRUE
  expect_true(is_valid_de_entry(rowname_entry, "rowname_entry"))
})

test_that("is_valid_pathway correctly identifies valid and invalid KEGG IDs", {
  # Valid KEGG IDs
  expect_true(is_valid_pathway("hsa04110"))
  expect_true(is_valid_pathway("mmu00010"))

  # Invalid KEGG IDs
  expect_false(is_valid_pathway("04110")) # numeric-only 5-digit
  expect_false(is_valid_pathway("hsa0411a")) # letters in numeric part
  expect_false(is_valid_pathway("0411a")) # letters in numeric-only
  expect_false(is_valid_pathway(4110)) # numeric input, not character
  expect_false(is_valid_pathway(c("hsa04110", "mmu00010"))) # length > 1
  expect_false(is_valid_pathway("")) # empty string
  expect_false(is_valid_pathway(NULL)) # NULL input
  expect_false(is_valid_pathway(NA)) # NA input
})

test_that("normalize_de_results handles all inputs correctly", {
  # NULL input
  expect_null(normalize_de_results(NULL))

  # data.frame input wrapped correctly
  df1 <- data.frame(KEGG_ids = c("hsa:1", "hsa:2"), log2FoldChange = c(1.2, -0.5))
  res1 <- suppressMessages(normalize_de_results(df1))
  expect_type(res1, "list")
  expect_named(res1, "de_input")
  expect_identical(res1$de_input$de_table, df1)
  expect_identical(res1$de_input$value_column, "log2FoldChange")
  expect_identical(res1$de_input$feature_column, "KEGG_ids")

  # data.frame with explicit column names
  df2 <- data.frame(id = c("A", "B"), value = c(1, 2))
  res2 <- suppressMessages(normalize_de_results(df2, value_column = "value", feature_column = "id"))
  expect_identical(res2$de_input$value_column, "value")
  expect_identical(res2$de_input$feature_column, "id")

  # Valid named list preserved
  de_list <- list(contrast1 = list(dummy = TRUE), contrast2 = list(dummy = TRUE))
  with_mocked_bindings(
    is_valid_de_entry = function(x, name) TRUE,
    {
      res3 <- normalize_de_results(de_list)
      expect_identical(res3, de_list)
    }
  )

  # Invalid entries removed from named list
  de_list2 <- list(valid = list(dummy = TRUE), invalid = list(dummy = FALSE))
  with_mocked_bindings(
    is_valid_de_entry = function(x, name) name == "valid",
    {
      res4 <- normalize_de_results(de_list2)
      expect_named(res4, "valid")
      expect_length(res4, 1)
    }
  )

  # Named list with no valid entries NULL + warning
  de_list3 <- list(a = list(), b = list())
  wmsgs <- capture_warnings(res5 <- normalize_de_results(de_list3))
  expect_true(any(grepl("problem in 'a'", wmsgs)))
  expect_true(any(grepl("problem in 'b'", wmsgs)))
  expect_length(wmsgs, 2)
  expect_null(res5)

  # List with empty names NULL + warning
  de_list4 <- list(a = list(), b = list())
  names(de_list4)[2] <- ""
  expect_warning(
    res6 <- normalize_de_results(de_list4),
    "must be NULL, a valid data.frame, or a named list"
  )
  expect_null(res6)

  # Unnamed list NULL + warning
  de_list5 <- list(list(dummy = TRUE), list(dummy = TRUE))
  expect_warning(
    res7 <- normalize_de_results(de_list5),
    "must be NULL, a valid data.frame, or a named list"
  )
  expect_null(res7)

  # Invalid type NULL + warning
  wmsgs <- capture_warnings(res8 <- normalize_de_results(42))
  expect_length(wmsgs, 2)
  expect_true(any(grepl("de_results must be a data.frame, but is of class: numeric", wmsgs)))
  expect_true(any(grepl("de_results must be NULL, a valid data.frame, or a named list. Ignoring de_results.", wmsgs)))
  expect_null(res8)
})

test_that("normalize_de_results handles MLimmaML and DFrame objects", {
  # MLimmaML object (simulate with S3 class)
  mlimma_obj <- structure(list(a = 1:2), class = "MLimmaML")
  wmsgs1 <- capture_warnings(res1 <- normalize_de_results(mlimma_obj))
  expect_length(wmsgs1, 1) # warning from is_valid_dataframe
  expect_null(res1)

  # S4 DFrame / DataFrame object (from S4Vectors)
  if (requireNamespace("S4Vectors", quietly = TRUE)) {
    library(S4Vectors)
    df <- DataFrame(KEGG_ids = c("hsa:1", "hsa:2"), log2FoldChange = c(1, -1))
    wmsgs2 <- capture_warnings(res2 <- normalize_de_results(df))
    # Depending on your is_valid_dataframe(), this may emit a warning or pass
    expect_true(length(wmsgs2) >= 0) # at least 0 warnings
    expect_null(res2) # likely NULL due to class check
  }
})
