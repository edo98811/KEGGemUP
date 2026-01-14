test_that("combine_results_in_dataframe correctly merges DE results across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    result <- combine_results_in_dataframe(de_list)

    # Structure checks
    expect_true(is.data.frame(result), info = test_name)
    expect_equal(colnames(result), c("ids_for_mapping", "de_value", "source"),
      info = paste0(test_name, " unexpected column names")
    )

    # Source tracking
    expect_setequal(unique(result$source), names(de_list))

    # IDs cleanup
    expect_false(any(grepl("^hsa:", result$ids_for_mapping)), info = test_name)
    expect_false(any(grepl("^cpd:", result$ids_for_mapping)), info = test_name)
    expect_false(any(grepl("^path:", result$ids_for_mapping)), info = test_name)

    # Missing value checks
    expect_false(any(is.na(result$ids_for_mapping)), info = test_name)
    expect_false(any(is.na(result$de_value)), info = test_name)

    # Row count check
    expected_nrows <- sum(vapply(de_list, function(x) nrow(x$de_table), numeric(1)))
    expect_equal(nrow(result), expected_nrows,
      info = paste0(test_name, " unexpected number of merged rows")
    )
  })
})

test_that("add_results_nodes correctly maps DE results onto nodes_df across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    results_combined <- combine_results_in_dataframe(de_list)

    if (test_name %in% throw_warning) {
      expect_warning(mapped_nodes <- add_results_nodes(expected_nodes, results_combined))
    } else {
      mapped_nodes <- add_results_nodes(expected_nodes, results_combined)
    }

    # Structure checks
    expect_true(is.data.frame(mapped_nodes), info = test_name)
    expect_true(all(c("id", "de_value", "color", "source", "text") %in% colnames(mapped_nodes)),
      info = paste0(test_name, " missing expected columns")
    )
    expect_equal(nrow(mapped_nodes), nrow(expected_nodes),
      info = paste0(test_name, " wrong number of rows")
    )

    # Content checks
    expect_false(all(is.na(mapped_nodes$de_value)), info = test_name)
    expect_false(all(is.na(mapped_nodes$source)), info = test_name)
    expect_false(all(is.na(mapped_nodes$text)), info = test_name)
  })
})

test_that("add_colors_to_nodes assigns colors based on de_value across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    results_combined <- combine_results_in_dataframe(de_list)

    if (test_name %in% throw_warning) {
      expect_warning(mapped_nodes <- add_results_nodes(expected_nodes, results_combined))
    } else {
      mapped_nodes <- add_results_nodes(expected_nodes, results_combined)
    }

    colored_nodes <- add_colors_to_nodes(mapped_nodes)

    # Structural checks
    expect_true(is.data.frame(colored_nodes), info = test_name)
    expect_true("color" %in% colnames(colored_nodes),
      info = paste0(test_name, " missing 'color' column")
    )

    # Color validity checks
    non_na_colors <- colored_nodes$color[!is.na(colored_nodes$color)]
    expect_false(all(is.na(colored_nodes$color)), info = test_name)
    expect_true(all(grepl("^#([A-Fa-f0-9]{6})$", non_na_colors)),
      info = paste0(test_name, " invalid color hex codes detected")
    )
  })
})
