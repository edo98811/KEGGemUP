library(mockery)

test_that("expand_keggs handles multiple KEGG IDs correctly", {
  input <- data.frame(
    name = c(1, 2, 3),
    KEGG = c("hsa:1234;hsa:5678", "cpd:C00022", "ko:K00001;ko:K00002;ko:K00003"),
    stringsAsFactors = FALSE
  )
  expected_output <- data.frame(name = c(1, 1, 2, 3, 3, 3), KEGG = c("1234", "5678", "C00022", "K00001", "K00002", "K00003"))
  actual_output <- expand_keggs(input)
  expect_equal(actual_output, expected_output)
})

test_that("remove_kegg_prefix_str removes prefixes and handles multiple IDs", {
  input <- c("hsa:1234 hsa:5678", "cpd:C00022", "ko:K00001 ko:K00002 ko:K00003")
  expected_output <- c("1234;5678", "C00022", "K00001;K00002;K00003")
  actual_output <- vapply(input, remove_kegg_prefix_str, FUN.VALUE = character(1), USE.NAMES = FALSE)
  expect_equal(actual_output, expected_output)
})

test_that("parse_kgml_edges load relationsps correctly", {
  edges_df <- suppressMessages(parse_kgml_relations(kgml_path))

  expect_equal(edges_df, expected_edges)
})

test_that("parse_kgml_entries loads empty edges  correctly", {
  expect_warning(edges_df <- suppressMessages(parse_kgml_relations(kgml_path_empty)), "No relations found in KGML file.")

  expect_true(nrow(edges_df) == 0)
  expect_true(inherits(edges_df, "data.frame"))
  expect_equal(colnames(edges_df), c("from", "to", "type", "subtype", "rel_value"))
})

test_that("parse_kgml_entries load nodes correctly", {
  nodes_df <- suppressMessages(parse_kgml_entries(kgml_path))

  expect_equal(nodes_df, expected_nodes)
})

test_that("combine_results_in_dataframe correctly merges DE results across all test lists", {
  purrr::imap(all_de_test_lists, function(de_list, test_name) {
    # Run the function for this test list
    result <- combine_results_in_dataframe(de_list)

    # Structure checks
    expect_true(is.data.frame(result), info = test_name)
    expect_equal(colnames(result), c("KEGG", "plot_value", "source"),
      info = paste0(test_name, "unexpected column names")
    )

    # Source tracking
    expect_setequal(unique(result$source), names(de_list))

    # KEGG ID cleanup
    expect_false(any(grepl("^hsa:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^cpd:", result$KEGG)), info = test_name)
    expect_false(any(grepl("^path:", result$KEGG)), info = test_name)

    # Missing value checks
    expect_false(any(is.na(result$KEGG)), info = test_name)
    expect_false(any(is.na(result$plot_value)), info = test_name)

    # Row count check
    expected_nrows <- sum(vapply(de_list, function(x) nrow(x$de_table), numeric(1)))
    expect_equal(nrow(result), expected_nrows,
      info = paste0(test_name, "unexpected number of merged rows")
    )
  })
})

test_that("add_compound_names caches and assigns glycan and compounds names", {
  bfc <- BiocFileCache(tempfile(), ask = FALSE)

  with_mocked_bindings(
    get_kegg_db = function(bfc_arg, db_type) {
      if (db_type == "compound") {
        return(data.frame(
          name = rownames(real_compounds),
          value = real_compounds[[1]],
          stringsAsFactors = FALSE
        ))
      } else if (db_type == "glycan") {
        return(data.frame(
          name = rownames(real_glycans),
          value = real_glycans[[1]],
          stringsAsFactors = FALSE
        ))
      } else {
        stop("Unexpected db_type")
      }
    },
    {
      res <- add_compound_names(nodes_df_basic, bfc)
    }
  )
  # res <- add_compound_names(nodes_df_basic, bfc)
  # gene node unchanged
  expect_equal(res$label[1], NA_character_)

  # known compound gets correct label

  expect_equal(res$label[2], gsub(";.*", "", as.character(real_compounds["C00001", 1])))
  expect_equal(res$label[4], gsub(";.*", "", as.character(real_glycans["G00001", 1])))

  # unknown compound keeps original ID
  expect_equal(res$label[3], "C99999")
})


test_that("download_kgml rejects invalid inputs", {
  expect_error(
    download_kgml("hsa00010"),
    "Either 'directory' or 'bfc' must be provided"
  )

  expect_error(
    download_kgml("hsa00010", bfc = 1),
    "BiocFileCache"
  )

  expect_error(
    download_kgml("hsa00010", directory = c("a", "b")),
    "single string"
  )
})

test_that("download_kgml works in directory mode", {
  tmpdir <- tempdir()
  expected_file <- file.path(tmpdir, "hsa00010.xml")

  result <- suppressMessages(download_kgml(
    pathway_id = "hsa00010",
    directory = tmpdir
  ))

  expect_equal(result, expected_file)
  expect_true(file.exists(result))
})

test_that("download_kgml works in cache mode", {
  fake_bfc <- BiocFileCache(tempdir(), ask = FALSE)

  result <- suppressMessages(download_kgml(
    pathway_id = "hsa00010",
    bfc = fake_bfc
  ))

  fake_path <- BiocFileCache::bfcrpath(fake_bfc, "https://rest.kegg.jp/get/hsa00010/kgml")
  expect_equal(result, fake_path)
})
