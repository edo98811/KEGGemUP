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
