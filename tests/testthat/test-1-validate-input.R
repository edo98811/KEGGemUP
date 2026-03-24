test_that("invalid inputs return NULL", {
  
  expect_null(standardize_de_results(NULL))
  
  expect_warning(
    res <- standardize_de_results(data.frame(a = 1)),
    "value_column and feature_column"
  )
  expect_null(res)
  
  expect_warning(
    res <- standardize_de_results(list(a = 1), value_column = "x"),
    "value_column and feature_column"
  )
  expect_null(res)
  
})


test_that("single dataframe input is standardized correctly", {
  
  df <- data.frame(
    gene = c("a", "b", "c"),
    logFC = c(1.2, -0.3, 0.8),
    stringsAsFactors = FALSE
  )
  
  res <- standardize_de_results(
    df,
    value_column = "logFC",
    feature_column = "gene"
  )
  
  expect_type(res, "list")
  expect_named(res, "de_input")
  
  entry <- res$de_input
  
  expect_true(is.data.frame(entry$de_table))
  expect_identical(entry$value_column, "logFC")
  expect_identical(entry$feature_column, "gene")
})


test_that("single dataframe with missing columns returns NULL", {
  
  df <- data.frame(
    gene = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    res <- standardize_de_results(
      df,
      value_column = "logFC",
      feature_column = "gene"
    ),
    "required columns"
  )
  
  expect_null(res)
  
})

test_that("single dataframe with wrong column types returns NULL", {
  
  df <- data.frame(
    gene = c("a", "b", "c"),
    logFC = c("1", "2", "3"),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    res <- standardize_de_results(
      df,
      value_column = "logFC",
      feature_column = "gene"
    )
  )
  
  expect_null(res)
  
})

test_that("valid list input is preserved", {
  
  df <- data.frame(
    gene = c("a", "b", "c"),
    logFC = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  
  input <- list(
    test1 = list(
      de_table = df,
      value_column = "logFC",
      feature_column = "gene"
    )
  )
  
  res <- standardize_de_results(input)
  
  expect_length(res, 1)
  expect_named(res, "test1")
  expect_true(is.data.frame(res$test1$de_table))
})


test_that("invalid list entries are removed", {
  
  df <- data.frame(
    gene = c("a", "b"),
    logFC = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  bad_df <- data.frame(
    gene = c("a", "b"),
    logFC = c("1", "2"),
    stringsAsFactors = FALSE
  )
  
  input <- list(
    valid = list(
      de_table = df,
      value_column = "logFC",
      feature_column = "gene"
    ),
    invalid_type = list(
      de_table = bad_df,
      value_column = "logFC",
      feature_column = "gene"
    ),
    missing_column = list(
      de_table = df,
      value_column = "missing",
      feature_column = "gene"
    )
  )
  
  res <- standardize_de_results(input)
  
  expect_length(res, 1)
  expect_named(res, "valid")
  
})

test_that("all invalid list entries return NULL", {
  
  df <- data.frame(
    gene = c("a", "b"),
    logFC = c("1", "2"),
    stringsAsFactors = FALSE
  )
  
  input <- list(
    bad = list(
      de_table = df,
      value_column = "logFC",
      feature_column = "gene"
    )
  )
  
  res <- standardize_de_results(input)
  
  expect_null(res)
  
})

test_that("list must be named", {
  
  df <- data.frame(
    gene = c("a", "b"),
    logFC = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  input <- list(
    list(
      de_table = df,
      value_column = "logFC",
      feature_column = "gene"
    )
  )
  
  expect_warning(
    res <- standardize_de_results(input)
  )
  
  expect_null(res)
  
})