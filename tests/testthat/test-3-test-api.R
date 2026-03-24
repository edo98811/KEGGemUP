# 
# test_that("download_kgml rejects invalid inputs", {
#   expect_error(
#     download_kgml("hsa00001", bfc = 1),
#     "BiocFileCache"
#   )
# 
#   expect_error(
#     download_kgml(
#       "hsa00001",
#       path = c("a", "b")),
#     "single string"
#   )
# })
# 
# test_that("download_kgml works in directory mode", {
#   tmpdir <- tempdir()
#   expected_file <- file.path(tmpdir, "hsa00010.xml")
# 
#   result <- suppressMessages(download_kgml(
#     pathway_id = "hsa00010",
#     path = tmpdir
#   ))
# 
#   expect_equal(result, expected_file)
#   expect_true(file.exists(result))
# })
# 
# test_that("download_kgml works in cache mode", {
#   result <- suppressMessages(download_kgml(
#     pathway_id = "hsa00010",
#     bfc = bfc
#   ))
# 
#   fake_path <- BiocFileCache::bfcrpath(bfc, "https://rest.kegg.jp/get/hsa00010/kgml")  # Changed to a real pathway ID
#   expect_equal(result, fake_path)
# })
# 
# test_that("get_kegg_db rejects invalid inputs", {
#   expect_error(
#     get_kegg_db("compound", bfc = 1),
#     "BiocFileCache"
#   )
# 
#   expect_error(
#     get_kegg_db("compound", path = c("a", "b")),
#     "single string"
#   )
# })
# 
# test_that("get_kegg_db works in directory mode", {
#   tmpdir <- tempdir()
#   expected_file <- file.path(tmpdir, "kegg_compound.tsv")
# 
#   result <- suppressMessages(
#     get_kegg_db(
#       db_name = "compound",
#       path = tmpdir
#     )
#   )
# 
#   expect_true(file.exists(expected_file))
#   expect_s3_class(result, "data.frame")
# })
