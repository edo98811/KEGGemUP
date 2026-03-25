test_that("combine_results_in_dataframe correctly merges DE results", {
  result <- KEGGemUP:::combine_results_in_dataframe(de_results_list)

  # Structure checks
  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("ids_for_mapping", "de_value", "de_source", "de_name"))

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

test_that("get_palette_range works correctly", {
  legend_input_1 <- c(0.5)
  legend_input_na <- c(0.5, NA, 0.8)
  legend_input_empty <- numeric(0)
  legend_input_n <- c(-2, -1, 0, 1, 2, 3)
  legend_input_inf <- c(-2, -1, Inf, 1, 2, 3)

  colors <- RColorBrewer::brewer.pal(11, "RdYlGn")
  palette_ramp <- colorRampPalette(colors)

  # Expected values
  suppressMessages({
    expect_true(is.numeric(KEGGemUP:::get_palette_range(legend_input_1, source_name = "test")))
    expect_true(is.numeric(KEGGemUP:::get_palette_range(legend_input_na, source_name = "test")))
    expect_true(is.na(KEGGemUP:::get_palette_range(legend_input_empty, source_name = "test")))
    expect_true(is.numeric(KEGGemUP:::get_palette_range(legend_input_n, source_name = "test")))

    expect_warning(res <- KEGGemUP:::get_palette_range(legend_input_inf, source_name = "test"))
    expect_true(is.na(res))
  })
})

test_that("add_results_nodes with malformed vertices_df", {
  # Missing ids_for_mapping column
  malformed_vertices_df <- data.frame(name = c("node2", "node2"))
  results_combined <- KEGGemUP:::combine_results_in_dataframe(de_results_list)
  expect_error(
    expect_warning(
      KEGGemUP:::add_results_nodes(malformed_vertices_df, results_combined)
    ),
    regexp = "add_results_nodes: Invalid input"
  )

  # ids_for_mapping column with non-character values
  malformed_vertices_df <- data.frame(name = c(1, 2), ids_for_mapping = c(1, 2))
  results_combined <- KEGGemUP:::combine_results_in_dataframe(de_results_list)

  expect_error(
    expect_warning(
      KEGGemUP:::add_results_nodes(malformed_vertices_df, results_combined)
    ),
    regexp = "add_results_nodes: Invalid input"
  )

  # ids_for_mapping column with NA
  malformed_vertices_df$ids_for_mapping <- c(NA, "test")
  expect_error(
    expect_warning(
      mapped_nodes <- KEGGemUP:::add_results_nodes(malformed_vertices_df, results_combined)
    ),
    regexp = "add_results_nodes: Invalid input"
  )
})

test_that("add_colors_to_nodes correctly maps DE results onto vertices_df", {
  # This is the test in case everything works correctly,
  # We prepare the vertices_df as it would be after parsing the KGML,
  # with the ids_for_mapping column ready for mapping.

  # Prepare vertices_df
  vertices_df <- kgml_steps$all_nodes
  indexes_to_map <- which(
    vertices_df$graphics_type != "line" & vertices_df$graphics_type != "group"
  )
  vertices_df$ids_for_mapping[indexes_to_map] <- vapply(
    vertices_df$KEGG[indexes_to_map],
    remove_kegg_prefix_str,
    character(1)
  )

  # Map results from the results list
  # We expect a warning because some nodes have multiple matching IDs
  results_combined <- KEGGemUP:::combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- KEGGemUP:::add_results_nodes(vertices_df, results_combined), "Some nodes had multiple matching IDs")

  # Expect dataframe with expected columns and same number of rows as input vertices_df
  expect_true(is.data.frame(mapped_nodes))
  expect_true(
    all(c("name", "de_value", "vertex.color", "de_source", "de_text", "de_name") %in% colnames(mapped_nodes)),
    info = "missing expected columns"
  )
  expect_equal(
    nrow(mapped_nodes),
    nrow(vertices_df),
    info = "wrong number of rows"
  )

  # Everything should have at least a non NA value
  expect_false(all(is.na(mapped_nodes$de_value)))
  expect_false(all(is.na(mapped_nodes$de_source)))
  expect_false(all(is.na(mapped_nodes$de_text)))
  expect_false(all(is.na(mapped_nodes$de_name)))
})

test_that("get_palette_colors works correctly and handles errors", {
  palette_valid <- "RdYlGn"
  palette_invalid <- "NotAValidPalette"
  palette_char_vector <- c("#FF0000", "#00FF00", "#0000FF")
  palette_null <- NULL

  # Test with valid RColorBrewer palette
  expect_true(is.character(KEGGemUP:::get_palette_colors(palette_valid)))

  # Test with invalid RColorBrewer palette
  expect_warning(
    KEGGemUP:::get_palette_colors(palette_invalid),
    "Palette 'NotAValidPalette' is not a valid RColorBrewer palette."
  )

  # Test with character vector of colors
  expect_equal(KEGGemUP:::get_palette_colors(palette_char_vector), palette_char_vector)
})


test_that("validate_palettes handles different scenarios", {
  sources <- c("A", "B")

  # Palettes_list wrong type (data.frame)
  palettes_list <- data.frame(x = 1:2)
  expect_error(
    result <- KEGGemUP:::validate_palettes(sources, palettes_list),
    "`palettes_list` must be a character vector"
  )

  # Palettes_list wrong names
  palettes_list <- c(X = "Spectral", Y = "Viridis")
  expect_warning(
    result <- KEGGemUP:::validate_palettes(sources, palettes_list),
    "Some sources do not have specified palettes"
  )
  expect_equal(
    result,
    setNames(
      list(
        rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")),
        rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))
      ),
      sources
    )
  )

  # Palette with global invalid palette
  warn <- capture_warnings(
    result <- KEGGemUP:::validate_palettes(sources, palettes_list, palette = "test"),
  )
  expect_equal(length(warn), 2)
  expect_equal(
    result,
    setNames(
      list(
        rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")),
        rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))
      ),
      sources
    )
  )
  
  # Palette with global palette and palettes list with one invalid
  warn <- capture_warnings(
    result <- KEGGemUP:::validate_palettes(sources, palettes_list <- c(X = "Spectral", A = "Viridis"), palette = "test")
  )
  expect_equal(length(warn), 2)

  # Correct list with proper names
  palettes_list <- c(A = "Spectral", B = "RdPu")
  result <- KEGGemUP:::validate_palettes(sources, palettes_list)
  expect_equal(
    result,
    setNames(
      list(
        rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")),
        rev(RColorBrewer::brewer.pal(n = 7, name = "RdPu"))
      ),
      sources
    )
  )

  # Palettes_list empty with custom default
  result <- KEGGemUP:::validate_palettes(sources, palette = "RdPu")
  expect_equal(
    result,
    setNames(
      list(
        rev(RColorBrewer::brewer.pal(n = 7, name = "RdPu")),
        rev(RColorBrewer::brewer.pal(n = 7, name = "RdPu"))
      ),
      sources
    )
  )

  # Test with invalid palette
  expect_warning(
    a <- KEGGemUP:::validate_palettes(c("A"), palette = "invalid"),
    "Palette 'invalid'"
  )
})

test_that("validate_palette_limits handles different scenarios", {
  sources <- c("A", "B")

  # Palettes_limits_list wrong type
  palettes_limits_list <- "Test"
  expect_error(
    result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list, palette_limit = NULL, default_palette_limit = FALSE),
    "`palettes_limits_list` must be a numeric vector"
  )

  # Palettes_limits_list wrong names
  palettes_limits_list <- c(X = 0.5, Y = 1)
  expect_warning(
    result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list, palette_limit = NULL, default_palette_limit = FALSE),
    "Some sources do not have specified palette limits"
  )
  expect_equal(result, setNames(rep(FALSE, 2), sources))

  # Correct list with proper names
  palettes_limits_list <- c(A = 0.5, B = 1)
  result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list = palettes_limits_list, palette_limit = NULL, default_palette_limit = FALSE)
  expect_equal(result, palettes_limits_list)

  # Palettes_limits_list not set returns default
  result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list = c(NA_real_), palette_limit = NULL, default_palette_limit = 1)
  expect_equal(result, setNames(rep(1, 2), sources))

  # Palettes_limits_list empty with custom default
  result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list = c(NA_real_), palette_limit = NULL, default_palette_limit = FALSE)
  expect_equal(result, setNames(rep(FALSE, 2), sources))

  # Test with palette limit and no list
  result <- KEGGemUP:::validate_palette_limits(sources, palettes_limits_list = c(NA_real_), palette_limit = 2, default_palette_limit = FALSE)
  expect_equal(result, setNames(rep(2, 2), sources))
})

test_that("add_colors_to_nodes assigns colors based on de_value", {
  # Same as before, (would it make sense to put this in a function?)

  # Prepare vertices_df
  vertices_df <- kgml_steps$all_nodes
  indexes_to_map <- which(
    vertices_df$graphics_type != "line" & vertices_df$graphics_type != "group"
  )
  vertices_df$ids_for_mapping[indexes_to_map] <- vapply(
    vertices_df$KEGG[indexes_to_map],
    KEGGemUP:::remove_kegg_prefix_str,
    character(1)
  )

  # Map results and add colors
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- KEGGemUP:::add_results_nodes(
    vertices_df,
    results_combined
  ), "Some nodes had multiple matching IDs")

  # Test with valid palette (get dataframe)
  suppressMessages(
    colored_nodes <- KEGGemUP:::add_colors_to_nodes(
      mapped_nodes,
      palette = "RdYlGn",
      palettes_limits_list = c(NA_real_),
      palette_limit = NULL,
      palettes_list = c(NA_character_)
    )$vertices_df
  )

  # Structural checks
  expect_true(is.data.frame(colored_nodes))
  expect_true("vertex.color" %in% colnames(colored_nodes))

  # Color assignment checks
  expect_true(any(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colored_nodes$vertex.color)))
  expect_true(any(colored_nodes$vertex.color != "white"))
})

test_that("add_colors_to_nodes handles errors in palette validation", {

  # Prepare vertices_df
  vertices_df <- kgml_steps$all_nodes
  indexes_to_map <- which(
    vertices_df$graphics_type != "line" & vertices_df$graphics_type != "group"
  )
  vertices_df$ids_for_mapping[indexes_to_map] <- vapply(
    vertices_df$KEGG[indexes_to_map],
    KEGGemUP:::remove_kegg_prefix_str,
    character(1)
  )

  # Map results and add colors with invalid palette
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- KEGGemUP:::add_results_nodes(
    vertices_df,
    results_combined
  ), "Some nodes had multiple matching IDs")
  
  suppressMessages(
    warns <- capture_warnings(
      colored_nodes <- KEGGemUP:::add_colors_to_nodes(
        mapped_nodes,
        palette = c("sasf", "invalid"),
        palettes_limits_list = c(NA_real_),
        palette_limit = NULL,
        palettes_list = c(NA_character_)
      )$vertices_df
    )
  )
  
  expect_equal(length(warns), 8)
  # Should still return a data frame with vertex.color column
  expect_true(is.data.frame(colored_nodes))
})