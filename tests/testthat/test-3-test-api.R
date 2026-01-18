test_that("add_node_labels caches and assigns glycan, compound, and gene names", {
  # Run the function
  res <- add_node_labels(nodes_df_basic, bfc = bfc)

  # gene node unchanged if not in database
  expect_equal(res$label[1], NA_character_)

  # known compound gets correct label
  expected_c00001 <- gsub(";.*", "", as.character(real_compounds["C00001", 1]))
  expect_equal(res$label[2], expected_c00001)

  # known glycan gets correct label
  expected_g00001 <- gsub(";.*", "", as.character(real_glycans["G00001", 1]))
  expect_equal(res$label[4], expected_g00001)

  # unknown compound keeps original ID
  expect_equal(res$label[3], "C99999")

  # gene node with known KEGG ID gets correct label
  expected_k00001 <- gsub(";.*", "", as.character(real_genes["K00001", 1]))
  expect_equal(res$label[5], expected_k00001)
})

test_that("add_reaction_labels assigns reaction labels and links correctly", {

  # Run the function
  res <- add_reaction_labels(nodes_df_basic, bfc = bfc)

  # reaction with known ID gets correct label
  expected_r00001 <- gsub(";.*", "", as.character(real_reactions["R00001", 1]))
  expect_equal(res$reaction_label[1], expected_r00001)

  # reaction with unknown ID keeps original
  expect_equal(res$reaction_label[2], "R99999")

  # reaction link created for known reaction
  expect_equal(
    res$reaction_link[1],
    "https://www.kegg.jp/dbget-bin/www_bget?R00001"
  )

  # reaction link is NA for unknown reaction
  expect_true(is.na(res$reaction_link[2]))
})


test_that("download_kgml rejects invalid inputs", {
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

  result <- suppressMessages(download_kgml(
    pathway_id = "hsa00010",
    bfc = fake_bfc
  ))

  fake_path <- BiocFileCache::bfcrpath(fake_bfc, "https://rest.kegg.jp/get/hsa00010/kgml")
  expect_equal(result, fake_path)
})

test_that("get_kegg_db rejects invalid inputs", {
  expect_error(
    get_kegg_db("compound", bfc = 1),
    "BiocFileCache"
  )

  expect_error(
    get_kegg_db("compound", directory = c("a", "b")),
    "single string"
  )
})

test_that("get_kegg_db works in directory mode", {
  tmpdir <- tempdir()
  expected_file <- file.path(tmpdir, "kegg_compound.tsv")

  result <- suppressMessages(
    get_kegg_db(
      db_name = "compound",
      directory = tmpdir
    )
  )

  expect_true(file.exists(expected_file))
  expect_s3_class(result, "data.frame")
})
