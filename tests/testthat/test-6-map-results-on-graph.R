test_that("combine_results_in_dataframe correctly merges DE results", {
  result <- combine_results_in_dataframe(de_results_list)

  # Structure checks
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("ids_for_mapping", "de_value", "de_source"))

  # Source tracking
  expect_setequal(unique(result$de_source), names(de_results_list))

  # Row count check
  expected_nrows <- sum(vapply(
    de_results_list,
    function(x) nrow(x$de_table),
    numeric(1)
  ))
  expect_equal(nrow(result), expected_nrows)
})

test_that("add_results_nodes correctly maps DE results onto nodes_df", {
  # Prepare nodes_df
  nodes_df <- kgml_steps$all_nodes
  indexes_to_map <- which(
    nodes_df$graphics_type != "line" & nodes_df$graphics_type != "group"
  )
  nodes_df$ids_for_mapping[indexes_to_map] <- vapply(
    nodes_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Map results
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- add_results_nodes(nodes_df, results_combined))

  # Structure checks
  expect_true(is.data.frame(mapped_nodes))
  expect_true(
    all(c("name", "de_value", "vertex.color", "de_source", "text") %in% colnames(mapped_nodes)),
    info = "missing expected columns"
  )
  expect_equal(
    nrow(mapped_nodes),
    nrow(nodes_df),
    info = "wrong number of rows"
  )

  # Content checks
  expect_false(all(is.na(mapped_nodes$de_value)))
  expect_false(all(is.na(mapped_nodes$de_source)))
  expect_false(all(is.na(mapped_nodes$text)))
})

test_that("add_colors_to_nodes assigns colors based on de_value", {

  # Prepare nodes_df
  nodes_df <- kgml_steps$all_nodes
  indexes_to_map <- which(
    nodes_df$graphics_type != "line" & nodes_df$graphics_type != "group"
  )
  nodes_df$ids_for_mapping[indexes_to_map] <- vapply(
    nodes_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Map results and add colors
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- add_results_nodes(
    nodes_df,
    results_combined
  ))
  colored_nodes <- add_colors_to_nodes(mapped_nodes)

  # Structural checks
  expect_true(is.data.frame(colored_nodes))
  expect_true("vertex.color" %in% colnames(colored_nodes))

  # Color assignment checks
  expect_true(any(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colored_nodes$vertex.color)))
  expect_true(any(colored_nodes$vertex.color != "white"))
})

test_that("add_results_nodes handles invalid ids_for_mapping column", {
  invalid_nodes_df <- kgml_steps$all_nodes
  results_combined <- combine_results_in_dataframe(de_results_list)

  # Test with missing ids_for_mapping column
  invalid_nodes_df$ids_for_mapping <- NA
  expect_warning(mapped_nodes <- add_results_nodes(invalid_nodes_df, results_combined))
  expect_true(is.data.frame(mapped_nodes))
  expect_equal(nrow(mapped_nodes), nrow(invalid_nodes_df))

  # Test with all NA ids_for_mapping
  invalid_nodes_df$ids_for_mapping <- ""
  expect_warning(mapped_nodes <- add_results_nodes(invalid_nodes_df, results_combined))
  expect_true(is.data.frame(mapped_nodes))
  expect_equal(nrow(mapped_nodes), nrow(invalid_nodes_df))

  invalid_nodes_df$ids_for_mapping <- NULL
  expect_error(
    add_results_nodes(invalid_nodes_df, results_combined),
    regexp = "Missing columns in nodes_df: ids_for_mapping"
  )
})
