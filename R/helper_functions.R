remove_kegg_prefix_str <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  separated_elements <- strsplit(kegg_ids, " ")
  ids <- lapply(separated_elements, function(x) sub("^[a-z]+:", "", x))
  ids <- paste(unlist(ids), collapse = ";")
  return(ids)
}

#' Convert KEGG IDs with prefixes to IDs without prefixes
#' @param kegg_ids Character vector of KEGG IDs with prefixes (e.g., "cpd:C00022", "mmu:1234", "ko:K00001 ko:K00002")
#' @return List of character vectors with KEGG IDs without prefixes
remove_kegg_prefix <- function(kegg_ids) {
  # Remove prefix (e.g., "cpd:", "mmu:", "ko:", "path:")
  ids <- vapply(kegg_ids, function(x) sub("^[a-z]+:", "", x), FUN.VALUE = character(1))
  return(ids)
}

#' Handle multiple KEGG IDs in a single string separated by ";"
#' @param kegg_df Data frame with at least two columns: 'name' and 'KEGG'
#' @return Data frame with a single column 'KEGG' containing individual KEGG IDs
expand_keggs <- function(kegg_df) {
  # Initialize empty vectors to store results
  ids_out <- c()
  kegg_out <- c()

  # Loop through each row of the data frame
  for (i in seq_len(nrow(kegg_df))) {
    # Split the KEGG string by ";"
    split_ids <- unlist(strsplit(kegg_df$KEGG[i], ";"))
    # Remove the prefix before ":" in each KEGG ID
    split_ids <- sub(".*:", "", split_ids)
    # Append the row ids and KEGG IDs
    ids_out <- c(ids_out, rep(kegg_df$name[i], length(split_ids)))
    kegg_out <- c(kegg_out, split_ids)
  }

  # Return the expanded data frame
  return(data.frame(name = ids_out, KEGG = kegg_out, stringsAsFactors = FALSE))
}

