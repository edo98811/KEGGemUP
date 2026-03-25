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

test_that("make_mapping_df works", {
  df1 <- data.frame(name = c("A","B"), ids_for_mapping = c("id1","id2"))
  res1 <- make_mapping_df(df1)
  expect_equal(nrow(res1), 2)
  expect_equal(res1$name, c("A","B"))
  expect_equal(res1$matched_id, c("id1","id2"))
  
  df2 <- data.frame(name = c("A","B"), ids_for_mapping = c("id1;id2","id3;id4;id5"))
  res2 <- make_mapping_df(df2)
  expect_equal(nrow(res2), 5)
  expect_equal(res2$name, c("A","A","B","B","B"))
  expect_equal(res2$matched_id, c("id1","id2","id3","id4","id5"))
  
  df3 <- data.frame(name = c("A","B","C"), ids_for_mapping = c("id1;id2","","id3"))
  res3 <- make_mapping_df(df3)
  expect_equal(nrow(res3), 3)
  expect_equal(res3$name, c("A","A","C"))
  expect_equal(res3$matched_id, c("id1","id2","id3"))
  
  df4 <- data.frame(name = c("A","B"), ids_for_mapping = c("",""))
  res4 <- make_mapping_df(df4)
  expect_equal(nrow(res4), NULL)
})