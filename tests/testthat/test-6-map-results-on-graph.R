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

test_that("get_palette_range works correctly") {
  legend_input_1 <- c(0.5)
  legend_input_na <- c(0.5, NA, 0.8)
  legend_input_empty <- numeric(0)
  legend_input_n <- c(-2, -1, 0, 1, 2, 3)
  legend_input_inf <- c(-2, -1, Inf, 1, 2, 3)

  # Expected values
  expect_true(is.numeric(create_legend_continous(legend_input_1, source = "test"))))
  expect_true(is.numeric(create_legend_continous(legend_input_na, source = "test"))))
  expect_true(is.na(create_legend_continous(legend_input_empty, source = "test"))))
  expect_true(is.numeric(create_legend_continous(legend_input_n, source = "test"))))
  expect_true(is.na(create_legend_continous(legend_input_inf, source = "test"))))
}

test_that("add_results_nodes with malformed vertices_df", {
  # Missing ids_for_mapping column
  malformed_vertices_df <- data.frame(name = c("node2", "node2"))
  results_combined <- combine_results_in_dataframe(de_results_list)

  expect_error(
    add_results_nodes(malformed_vertices_df, results_combined),
    regexp = "add_results_nodes: Invalid input"
  )
  # ids_for_mapping column with non-character values
  malformed_vertices_df <- data.frame(name = c(1, 2), ids_for_mapping = c(1, 2))
  results_combined <- combine_results_in_dataframe(de_results_list)

  expect_error(
    add_results_nodes(malformed_vertices_df, results_combined),
    regexp = "add_results_nodes: Invalid input"
  )

  # ids_for_mapping column with NA 
  malformed_vertices_df$ids_for_mapping <- c(NA, "test")
  expect_warning(
    mapped_nodes <- add_results_nodes(malformed_vertices_df, results_combined),
    regexp = "Some nodes had multiple matching IDs"
  )
  expect_true(is.data.frame(mapped_nodes))
  expect_equal(nrow(mapped_nodes), nrow(malformed_vertices_df))
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
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- add_results_nodes(vertices_df, results_combined),"Some nodes had multiple matching IDs")

  # Structure checks
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

  # Content checks
  expect_false(all(is.na(mapped_nodes$de_value)))
  expect_false(all(is.na(mapped_nodes$de_source)))
  expect_false(all(is.na(mapped_nodes$de_text)))
  expect_false(all(is.na(mapped_nodes$de_name)))
})

test_that("get_palette_colors "{
  palette_valid <- "RdYlGn"
  palette_invalid <- "NotAValidPalette"
  palette_char_vector <- c("#FF0000", "#00FF00", "#0000FF")
  palette_null <- NULL

  # Test with valid RColorBrewer palette
  expect_true(is.character(get_palette_colors(palette_valid)))

  # Test with invalid RColorBrewer palette
  expect_warning(get_palette_colors(palette_invalid),  "Palette 'NotAValidPalette' is not a valid RColorBrewer palette.")

  # Test with character vector of colors
  expect_equal(get_palette_colors(palette_char_vector), palette_char_vector)

  # Test with NULL palette
  expect_true(is.character(get_palette_colors(palette_null)))
})


test_that("validate_palettes handles different scenarios", {
  
  sources <- c("A", "B")
  
  # Palettes_list wrong type (data.frame)
  palettes_list <- data.frame(x = 1:2)
  expect_warning(
    result <- validate_palettes(sources, palettes_list),
    "Some sources do not have specified palettes"
  )
  expect_equal(result, setNames(rep("Spectral", 2), sources))
  
  # Palettes_list wrong names
  palettes_list <- list(X = "Spectral", Y = "Viridis")
  expect_warning(
    result <- validate_palettes(sources, palettes_list),
    "Some sources do not have specified palettes"
  )
  expect_equal(result, setNames(rep("Spectral", 2), sources))
  
  # Correct list with proper names
  palettes_list <- list(A = "Spectral", B = "Viridis")
  result <- validate_palettes(sources, palettes_list)
  expect_equal(result, palettes_list)
  
  # Palettes_list NULL returns default
  result <- validate_palettes(sources, palettes_list = NULL)
  expect_equal(result, setNames(rep("Spectral", 2), sources))
  
  # Palettes_list NULL with custom default
  result <- validate_palettes(sources, palettes_list = NULL, palette = "Viridis")
  expect_equal(result, setNames(rep("Viridis", 2), sources))
  
})

test_that("validate_palette_limits handles different scenarios", {
  
  sources <- c("A", "B")
  
  # Palettes_limits_list wrong type 
  palettes_limits_list <- "Test"
  expect_warning(
    result <- validate_palette_limits(sources, palettes_limits_list),
    "Some sources do not have specified palette limits"
  )
  expect_equal(result, setNames(rep(NULL, 2), sources))
  
  # Palettes_limits_list wrong names
  palettes_limits_list <- list(X = 0.5, Y = 1)
  expect_warning(
    result <- validate_palette_limits(sources, palettes_limits_list),
    "Some sources do not have specified palette limits"
  )
  expect_equal(result, setNames(rep(NULL, 2), sources))
  
  # Correct list with proper names
  palettes_limits_list <- list(A = 0.5, B = 1)
  result <- validate_palette_limits(sources, palettes_limits_list)
  expect_equal(result, palettes_limits_list)
  
  # Palettes_limits_list NULL returns default
  result <- validate_palette_limits(sources, palettes_limits_list = NULL, default_limit = 1)
  expect_equal(result, setNames(rep(1, 2), sources))
  
  # Palettes_limits_list NULL with custom default
  result <- validate_palette_limits(sources, palettes_limits_list = NULL, limit = 0.8)
  expect_equal(result, setNames(rep(0.8, 2), sources))
  
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
    remove_kegg_prefix_str,
    character(1)
  )

  # Map results and add colors
  results_combined <- combine_results_in_dataframe(de_results_list)
  expect_warning(mapped_nodes <- add_results_nodes(
    vertices_df,
    results_combined
  ), "Some nodes had multiple matching IDs")

  # Test with valid palette (get dataframe)
  colored_nodes <- add_colors_to_nodes(mapped_nodes,
  palette = "RdYlGn",)$vertices_df

  # Test with invalid palette
  expect_warning(
    add_colors_to_nodes(mapped_nodes, palette = c("invalid", "notacolor", "alsonotvalid")),
    "Failed to create color ramp for source 'invalid'"
  )

  # Structural checks
  expect_true(is.data.frame(colored_nodes))
  expect_true("vertex.color" %in% colnames(colored_nodes))

  # Color assignment checks
  expect_true(any(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colored_nodes$vertex.color)))
  expect_true(any(colored_nodes$vertex.color != "white"))
}
