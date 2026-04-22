# internal routines -------------------------------------------------------


#' Convert KEGG IDs with prefixes to a single string of IDs without prefixes
#'
#' @param kegg_ids Character string of KEGG IDs with prefixes
#' (e.g., 'cpd:C00022 mmu:1234 ko:K00001 ko:K00002')
#'
#' @return Character string of KEGG IDs without prefixes, separated by ';'
#'
#' @noRd
remove_kegg_prefix_str <- function(kegg_ids) {
  vapply(kegg_ids, function(id) {
    if (is.na(id)) {
      return(NA_character_)
    } # preserve NA
    elements <- strsplit(id, " ")[[1]] # split by space
    elements <- sub("^[a-z]+:", "", elements) # remove prefix
    paste(elements, collapse = ";") # collapse back to single string
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}

#' Convert KEGG IDs with prefixes to IDs without prefixes
#'
#' @param kegg_ids Character vector of KEGG IDs with prefixes
#' (e.g., 'cpd:C00022', 'mmu:1234', 'ko:K00001 ko:K00002')
#'
#' @return List of character vectors with KEGG IDs without prefixes
#'
#' @noRd
remove_kegg_prefix <- function(kegg_ids) {
  # Remove prefix (e.g., 'cpd:', 'mmu:', 'ko:', 'path:')
  ids <- vapply(kegg_ids, function(x) sub("^[a-z]+:", "", x), FUN.VALUE = character(1))
  return(ids)
}

#' Handle multiple KEGG IDs in a single string separated by ';'
#'
#' @param kegg_df Data frame with at least two columns: 'name' and 'KEGG'
#'
#' @return Data frame with a single column 'KEGG' containing individual KEGG IDs
#'
#' @noRd
expand_keggs <- function(kegg_df) {
  # Initialize empty vectors to store results
  ids_out <- c()
  kegg_out <- c()

  # Loop through each row of the data frame
  for (i in seq_len(nrow(kegg_df))) {
    # Split the KEGG string by ';'
    split_ids <- unlist(strsplit(kegg_df$KEGG[i], ";"))
    # Remove the prefix before ':' in each KEGG ID
    split_ids <- sub(".*:", "", split_ids)
    # Append the row ids and KEGG IDs
    ids_out <- c(ids_out, rep(kegg_df$name[i], length(split_ids)))
    kegg_out <- c(kegg_out, split_ids)
  }

  # Return the expanded data frame
  return(data.frame(name = ids_out, KEGG = kegg_out))
}


#' Create a mapping data frame
#'
#' Create a mapping data frame from vertices_df for subgraph extraction
#'
#' @param vertices_df Data frame of graph vertices with columns 'name' and
#' 'ids_for_mapping'
#'
#' @return A data frame with columns 'name' and 'matched_id'.
#' The 'matched_id' column contains individual KEGG IDs extracted from
#' 'ids_for_mapping', which may contain multiple IDs separated by ';'.
#'
#' @details Example of output:
#'   name matched_id
#' 1  NodeA     K00001
#' 2  NodeA     K00002
#' 3  NodeB     K00003
#'
#' @noRd
make_mapping_df <- function(vertices_df) {
  # Explode ids_for_mapping per node
  mapping <- do.call(
    rbind,
    lapply(seq_len(nrow(vertices_df)), function(i) {
      ids <- vertices_df$ids_for_mapping[i]
      if (ids == "") return(NULL)
      data.frame(
        name = vertices_df$name[i],
        matched_id = strsplit(
          ids, ";",
          fixed = TRUE
        )[[1]]
      )
    })
  )
  return(mapping)
}

