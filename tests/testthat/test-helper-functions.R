test_that("remove_kegg_prefix_str removes prefixes and preserves NA", {

  # Test vector with mixed IDs and NA
  input <- c("cpd:C00001 cpd:C00002", NA, "ko:K00001", "", "ec:1.1.1.1 ec:2.2.2.2")
  
  expected <- c(
    "C00001;C00002",  # prefixes removed, multiple IDs collapsed
    NA,               # NA preserved
    "K00001",         # prefix removed
    "",               # empty string stays empty
    "1.1.1.1;2.2.2.2" # multiple enzyme IDs collapsed
  )
  
  output <- remove_kegg_prefix_str(input)
  
  expect_equal(output, expected)
  
  # Confirm type of NA is still NA_character_
  expect_true(is.na(output[2]))
})