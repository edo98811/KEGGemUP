
test_that("add_results_nodes correctly maps DE results onto nodes_df across all test lists", {

  # Iterate through each DE test list
  purrr::imap(all_de_test_lists, function(de_list, test_name) {

    # Combine all DE results into one table
    results_combined <- combine_results_in_dataframe(de_list)

    # --- Capture all warnings ---
    warnings <- testthat::capture_warnings(
      mapped_nodes <- add_results_nodes(expected_nodes_cols, results_combined) # expect warnign expects only one not multiple!!
    )

    # --- Expected warning logic ---
    if (test_name %in% names(expected_warnings)) {
      exp_n <- expected_warnings[[test_name]]

      # Count and content validation
      expect_equal(length(warnings), exp_n,
                   info = paste0(test_name, " expected ", exp_n, " warnings, got, ", length(warnings)))
      expect_true(any(grepl("Multiple results", warnings)) || any(grepl("mapped", warnings)),
                  info = paste0(test_name, " expected warning about mapping issues"))
    } else {
      # Should produce no warnings
      expect_equal(length(warnings), 0,
                   info = paste0(test_name, " unexpected warnings found"))
    }

    # --- Structure checks ---
    expect_true(is.data.frame(mapped_nodes), info = test_name)

    expect_true(all(c("plot_value", "color", "source", "text") %in% colnames(mapped_nodes)),
      info = paste0(test_name, " missing expected columns")
    )
    expect_equal(nrow(mapped_nodes), nrow(expected_nodes_cols),
      info = paste0(test_name, " wrong number of rows")
    )

    # --- Content checks ---
    expect_false(all(is.na(mapped_nodes$plot_value)), info = test_name)
    expect_false(all(is.na(mapped_nodes$source)), info = test_name)
    expect_false(all(is.na(mapped_nodes$text)), info = test_name)
  })
})


test_that("add_colors_to_nodes assigns colors based on values across all test lists", {
  # Load the node reference table once

  purrr::imap(all_de_test_lists, function(de_list, test_name) {

    # Prepare combined results and map onto nodes
    results_combined <- combine_results_in_dataframe(de_list)

    # --- Capture all warnings ---
    warnings <- testthat::capture_warnings(
      mapped_nodes <- add_results_nodes(expected_nodes_cols, results_combined)
    )

    # --- Expected warning logic ---
    if (test_name %in% names(expected_warnings)) {
      exp_n <- expected_warnings[[test_name]]

      # Count and content validation
      expect_equal(length(warnings), exp_n,
        info = paste0(test_name, " expected ", exp_n, " warnings, got, ", length(warnings))
      )
      expect_true(any(grepl("Multiple results", warnings)) || any(grepl("mapped", warnings)),
        info = paste0(test_name, " expected warning about mapping issues")
      )
    } else {
      # Should produce no warnings
      expect_equal(length(warnings), 0,
        info = paste0(test_name, " unexpected warnings found")
      )
    }

    # Apply color assignment
    colored_nodes <- add_colors_to_nodes(mapped_nodes)

    # --- Structural checks ---
    expect_true(is.data.frame(colored_nodes), info = test_name)
    expect_true("color" %in% colnames(colored_nodes),
      info = paste0(test_name, "missing 'color' column")
    )

    # --- Color validity checks ---
    non_na_colors <- colored_nodes$color[!is.na(colored_nodes$color)]
    expect_false(all(is.na(colored_nodes$color)), info = test_name)
    expect_true(all(grepl("^#([A-Fa-f0-9]{6})$", non_na_colors)),
      info = paste0(test_name, "invalid color hex codes detected")
    )

    # --- Logical consistency ---
    expect_equal(!is.na(colored_nodes$plot_value), colored_nodes$color != "#FFFFFF")
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

  styled <- style_nodes(expected_nodes_cols)

  # expect_true(all(styled$fixed))
  expect_true(is.numeric(styled$widthConstraint))
  expect_true(is.numeric(styled$heightConstraint))
  expect_true(is.character(styled$label))
  expect_true(any(!is.na(styled$label)))

  expect_true(is.logical(styled$fixed))
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
    nodes_df_expected <- add_columns_nodes_df(nodes_df_expected)
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

test_that("add groups add groups correctly", {
  styled <- add_group(expected_nodes_cols)

  expect_equal(unique(styled$group), c("group_5", NA))
  expect_equal(sum(styled$group[!is.na(styled$group)] == "group_5"), 3)

})