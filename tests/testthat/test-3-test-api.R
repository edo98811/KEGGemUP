
test_that("download_kgml rejects invalid inputs", {
  expect_error(
    download_kgml("hsa00001", bfc = 1),
    "BiocFileCache"
  )

  expect_error(
    download_kgml(
      "hsa00001",
      path = c("a", "b")),
    "single string"
  )
})

test_that("download_kgml works in directory mode", {
  tmpdir <- tempdir()
  expected_file <- file.path(tmpdir, "hsa00010.xml")

  result <- suppressMessages(download_kgml(
    pathway_id = "hsa00010",
    path = tmpdir
  ))

  expect_equal(result, expected_file)
  expect_true(file.exists(result))
})

test_that("download_kgml works in cache mode", {
  result <- suppressMessages(download_kgml(
    pathway_id = "hsa00010",
    bfc = bfc
  ))

  fake_path <- BiocFileCache::bfcrpath(bfc, "https://rest.kegg.jp/get/hsa00010/kgml")  # Changed to a real pathway ID
  expect_equal(result, fake_path)
})

test_that("get_kegg_db rejects invalid inputs", {
  expect_error(
    get_kegg_db("compound", bfc = 1),
    "BiocFileCache"
  )

  expect_error(
    get_kegg_db("compound", path = c("a", "b")),
    "single string"
  )
})

test_that("get_kegg_db works in directory mode", {
  tmpdir <- tempdir()
  expected_file <- file.path(tmpdir, "kegg_compound.tsv")

  result <- suppressMessages(
    get_kegg_db(
      db_name = "compound",
      path = tmpdir
    )
  )

  expect_true(file.exists(expected_file))
  expect_s3_class(result, "data.frame")
})

test_that("select_cache_or_path works correctly", {
  
  # both provided -> error
  expect_error(
    select_cache_or_path(bfc = BiocFileCache::BiocFileCache(tempdir()), path = tempdir()),
    "Provide either 'bfc' OR 'path'"
  )
  
  # valid BiocFileCache -> cache mode
  bfc <- BiocFileCache::BiocFileCache(tempdir(), ask = FALSE)
  expect_equal(select_cache_or_path(bfc = bfc, path = NULL), "cache")
  
  # invalid bfc object -> error
  expect_error(
    select_cache_or_path(bfc = "not_a_cache", path = NULL),
    "'bfc' must be a valid BiocFileCache object"
  )
  
  # valid path -> dir mode
  p <- file.path(tempdir(), "kegg_test_dir")
  if (dir.exists(p)) unlink(p, recursive = TRUE)
  expect_equal(select_cache_or_path(bfc = NULL, path = p), "dir")
  expect_true(dir.exists(p))

  # null path -> none mode
  expect_equal(select_cache_or_path(bfc = NULL, path = NULL), "none")
  
  # invalid path type -> error
  expect_error(
    select_cache_or_path(bfc = NULL, path = c("a", "b")),
    "'path' must be a single string"
  )
  
})