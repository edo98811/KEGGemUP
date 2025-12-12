test_that("is_valid_de_entry works as expected", {
  # ---- Valid input ----
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

  # ---- Missing required list elements ----
  invalid_missing <- list(de_table = de_table)
  expect_warning(expect_false(is_valid_de_entry(invalid_missing, "invalid_missing")))

  # ---- de_table not a data frame ----
  invalid_table <- list(
    de_table = matrix(1:4, ncol = 2),
    value_column = "logFC",
    feature_column = "gene"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_table, "invalid_table")))

  # ---- value_column not present ----
  invalid_value_column <- list(
    de_table = de_table,
    value_column = "not_here",
    feature_column = "gene"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_value_column, "invalid_value_column")))

  # ---- feature_column not present or not rownames ----
  invalid_feature_column <- list(
    de_table = de_table,
    value_column = "logFC",
    feature_column = "not_here"
  )
  expect_warning(expect_false(is_valid_de_entry(invalid_feature_column, "invalid_feature_column")))

  # ---- feature_column = 'rownames' case ----
  rownames(de_table) <- de_table$gene
  rowname_entry <- list(
    de_table = de_table,
    value_column = "logFC",
    feature_column = "rownames"
  )

  # This checks logic– if 'rownames' allowed, expect TRUE
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
