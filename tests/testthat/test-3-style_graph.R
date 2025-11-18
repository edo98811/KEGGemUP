# test_that("add_results_nodes correctly maps DE results onto nodes_df", {
#   # Capture warnings to test multiple-match behavior later
#   nodes_df <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
#   results_combined <- combine_results_in_dataframe(de_results_list)
#   nodes_df$KEGG <- vapply(nodes_df$kegg_name, remove_kegg_prefix_str, FUN.VALUE = character(1))

#   # Run function
#   mapped_nodes <- add_results_nodes(nodes_df, results_combined)
#   expect_true(is.data.frame(mapped_nodes))
#  
#   # Check structure
#   expect_true(is.data.frame(mapped_nodes))
#   expect_true(all(c("value", "color", "source", "text") %in% colnames(mapped_nodes)))
#   expect_equal(nrow(mapped_nodes), nrow(nodes_df))
#   expect_false(all(is.na(mapped_nodes$value)))
#   expect_false(all(is.na(mapped_nodes$source)))
#   expect_false(all(is.na(mapped_nodes$text)))

#   results_combined_double <- combine_results_in_dataframe(de_results_list_with_repetition)
#   expect_warning(
#     add_results_nodes(nodes_df, results_combined_double),
#     "Multiple results mapped to node"
#   )
# })
test_that("add_results_nodes correctly maps DE results onto nodes_df across all test lists", {
  # Load node reference table
  nodes_df <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
  nodes_df$KEGG <- vapply(nodes_df$kegg_name, remove_kegg_prefix_str, FUN.VALUE = character(1))

  # Iterate through each DE test list
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    message(glue::glue("🔍 Testing mapping for scenario: {test_name}"))

    # Combine all DE results into one table
    results_combined <- combine_results_in_dataframe(de_list)

    # --- Capture all warnings ---
    warnings <- testthat::capture_warnings(
      mapped_nodes <- add_results_nodes(nodes_df, results_combined) # expect warnign expects only one not multiple!!
    )

    # --- Expected warning logic ---
    if (test_name %in% names(expected_warnings)) {
      exp_n <- expected_warnings[[test_name]]

      # Count and content validation
      expect_equal(length(warnings), exp_n,
                   info = glue::glue("{test_name}: expected {exp_n} warnings, got {length(warnings)}"))
      expect_true(any(grepl("Multiple results", warnings)) || any(grepl("mapped", warnings)),
                  info = glue::glue("{test_name}: expected warning about mapping issues"))
    } else {
      # Should produce no warnings
      expect_equal(length(warnings), 0,
                   info = glue::glue("{test_name}: unexpected warnings found"))
    }

    # --- Structure checks ---
    expect_true(is.data.frame(mapped_nodes), info = test_name)
    expect_true(all(c("value", "color", "source", "text") %in% colnames(mapped_nodes)),
      info = glue::glue("{test_name}: missing expected columns")
    )
    expect_equal(nrow(mapped_nodes), nrow(nodes_df),
      info = glue::glue("{test_name}: wrong number of rows")
    )

    # --- Content checks ---
    expect_false(all(is.na(mapped_nodes$value)), info = test_name)
    expect_false(all(is.na(mapped_nodes$source)), info = test_name)
    expect_false(all(is.na(mapped_nodes$text)), info = test_name)
  })
})


# test_that("add_color_to_nodes assigns colors based on values", {
#   nodes_df <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
#   results_combined <- combine_results_in_dataframe(de_results_list)
#   nodes_df$KEGG <- vapply(nodes_df$kegg_name, remove_kegg_prefix_str, FUN.VALUE = character(1))

#   nodes_df <- add_results_nodes(nodes_df, results_combined)

#   # Run function
#   nodes_df <- add_colors_to_nodes(nodes_df)

#   expect_true(is.data.frame(nodes_df))
#   expect_false(all(is.na(nodes_df$color)))
#   # Check non-NA colors are valid hex (#RRGGBB)
#   expect_true(all(grepl("^#([A-Fa-f0-9]{6})$", nodes_df$color[!is.na(nodes_df$color)])))
#   expect_equal(!is.na(nodes_df$value), !is.na(nodes_df$color))
# })

test_that("add_colors_to_nodes assigns colors based on values across all test lists", {
  # Load the node reference table once
  base_nodes_df <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
  base_nodes_df$KEGG <- vapply(base_nodes_df$kegg_name, remove_kegg_prefix_str, FUN.VALUE = character(1))

  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    message(glue::glue("🎨 Testing add_colors_to_nodes() for scenario: {test_name}"))

    # Prepare combined results and map onto nodes
    results_combined <- combine_results_in_dataframe(de_list)

    # --- Capture all warnings ---
    warnings <- testthat::capture_warnings(
      mapped_nodes <- add_results_nodes(base_nodes_df, results_combined)
    )

    # --- Expected warning logic ---
    if (test_name %in% names(expected_warnings)) {
      exp_n <- expected_warnings[[test_name]]

      # Count and content validation
      expect_equal(length(warnings), exp_n,
        info = glue::glue("{test_name}: expected {exp_n} warnings, got {length(warnings)}")
      )
      expect_true(any(grepl("Multiple results", warnings)) || any(grepl("mapped", warnings)),
        info = glue::glue("{test_name}: expected warning about mapping issues")
      )
    } else {
      # Should produce no warnings
      expect_equal(length(warnings), 0,
        info = glue::glue("{test_name}: unexpected warnings found")
      )
    }

    # Apply color assignment
    colored_nodes <- add_colors_to_nodes(mapped_nodes)

    # --- Structural checks ---
    expect_true(is.data.frame(colored_nodes), info = test_name)
    expect_true("color" %in% colnames(colored_nodes),
      info = glue::glue("{test_name}: missing 'color' column")
    )

    # --- Color validity checks ---
    non_na_colors <- colored_nodes$color[!is.na(colored_nodes$color)]
    expect_false(all(is.na(colored_nodes$color)), info = test_name)
    expect_true(all(grepl("^#([A-Fa-f0-9]{6})$", non_na_colors)),
      info = glue::glue("{test_name}: invalid color hex codes detected")
    )

    # --- Logical consistency ---
    expect_equal(!is.na(colored_nodes$value), colored_nodes$color != "#FFFFFF")
  })
})


test_that("known subtypes are styled correctly", {
  edges <- data.frame(
    subtype = c("activation", "inhibition", "phosphorylation"),
    id = 1:3
  )

  styled <- style_edges(edges)

  expect_equal(styled$dashes, c(FALSE, FALSE, FALSE))
  expect_equal(styled$arrows, c("to", "tee", "to"))
  expect_equal(styled$label, c("", "", "+p"))
})

test_that("unknown or NA subtypes default to others_unknown", {
  edges <- data.frame(
    subtype = c("nonsense", NA)
  )

  styled <- style_edges(edges)

  expect_true(all(styled$subtype == "others_unknown"))
  expect_true(all(styled$color == "black"))
  expect_true(all(styled$dashes))
  expect_true(all(styled$arrows == "to"))
  expect_true(all(styled$label == "?"))
})

test_that("slashes in subtype names are replaced with underscores", {
  edges <- data.frame(
    subtype = c("phosphorylation/dephosphorylation", "unknown/type")
  )

  styled <- style_edges(edges)
  expect_false(all(grepl("/", styled$subtype)))
  expect_equal(styled$subtype, c("others_unknown", "others_unknown")) # not defined, defaults
})

test_that("edges_df is empty or has one row", {
  edges_empty <- data.frame(subtype = character(0))
  styled_empty <- style_edges(edges_empty)
  expect_equal(nrow(styled_empty), 0)

  edges_one <- data.frame(subtype = c("activation"))
  styled_one <- style_edges(edges_one)
  expect_equal(nrow(styled_one), 1)
})

test_that("function returns same number of rows", {
  edges <- data.frame(subtype = c("activation", "repression", "state_change"))
  styled <- style_edges(edges)
  expect_equal(nrow(styled), 3)
})

test_that("adds visual styling columns correctly", {
  nodes <- data.frame(
    label = c("GeneA", "GeneB"),
    KEGG = c("hsa00010", "hsa00020"),
    type = c("compound", "simple"),
    height = 30,
    width = 30,
    x = 100,
    y = 200
  )

  styled <- style_nodes(nodes)

  expect_equal(styled$shape, c("dot", "box"))
  expect_true(all(styled$fixed))
  expect_true(all(c("x", "y", "widthConstraint", "heightConstraint") %in% names(styled)))
})

test_that("tooltip is formatted correctly", {
  nodes <- data.frame(
    label = "TP53",
    KEGG = "04115",
    type = "simple",
    kegg_name = "hsa04115",
    height = 30,
    width = 30,
    x = 100,
    y = 200
  )

  styled <- add_tooltip(nodes)
  expect_true("title" %in% names(styled))
})

test_that("make igraph works", {
    nodes_df_expected <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
    edges_df_expected <- tibble::as_tibble(read.csv(edges_df_path, sep = ";", colClasses = "character"))

    title = "Test Pathway"

    g_expected <- make_igraph_graph(nodes_df_expected, edges_df_expected, title)
    expect_true(igraph::is_igraph(g_expected))

    expect_warning(
      make_igraph_graph(nodes_df_expected, data.frame(from=character(0), to=character(0)), title),
      "No edges in graph."
    )

    expect_warning(
      make_igraph_graph(nodes_df_expected, NULL, title),
      "No edges in graph."
    )
})

test_that("make visnetwork graph works", {
    nodes_df_expected <- tibble::as_tibble(read.csv(nodes_df_path, sep = ";", colClasses = "character"))
    edges_df_expected <- tibble::as_tibble(read.csv(edges_df_path, sep = ";", colClasses = "character"))

    title = "Test Pathway"

    g_expected <- make_vis_graph(nodes_df_expected, edges_df_expected, title)
    expect_true(inherits(g_expected, "visNetwork"))

    expect_warning(
      make_vis_graph(nodes_df_expected, data.frame(from=character(0), to=character(0)), title),
      "No edges in graph."
    )
    expect_warning(
      make_vis_graph(nodes_df_expected, NULL, title),
      "No edges in graph."
    )
})
