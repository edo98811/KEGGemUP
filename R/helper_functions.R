#' Convert KEGG IDs with prefixes to a single string of IDs without prefixes
#' @param kegg_ids Character string of KEGG IDs with prefixes
#' (e.g., 'cpd:C00022 mmu:1234 ko:K00001 ko:K00002')
#' @return Character string of KEGG IDs without prefixes, separated by ';'
#' @noRd
remove_kegg_prefix_str <- function(kegg_ids) {
  # Remove prefix (e.g., 'cpd:', 'mmu:', 'ko:', 'path:')
  separated_elements <- strsplit(kegg_ids, " ")
  ids <- lapply(separated_elements, function(x) sub("^[a-z]+:", "", x))
  ids <- paste(unlist(ids), collapse = ";")
  return(ids)
}

#' Convert KEGG IDs with prefixes to IDs without prefixes
#' @param kegg_ids Character vector of KEGG IDs with prefixes
#' (e.g., 'cpd:C00022', 'mmu:1234', 'ko:K00001 ko:K00002')
#' @return List of character vectors with KEGG IDs without prefixes
#' @noRd
remove_kegg_prefix <- function(kegg_ids) {
  # Remove prefix (e.g., 'cpd:', 'mmu:', 'ko:', 'path:')
  ids <- vapply(kegg_ids, function(x) sub("^[a-z]+:", "", x), FUN.VALUE = character(1))
  return(ids)
}

#' Handle multiple KEGG IDs in a single string separated by ';'
#' @param kegg_df Data frame with at least two columns: 'name' and 'KEGG'
#' @return Data frame with a single column 'KEGG' containing individual KEGG IDs
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
  return(data.frame(name = ids_out, KEGG = kegg_out, stringsAsFactors = FALSE))
}

#' Return all cached KEGG and mapping files from BiocFileCache
#' @return A list containing data frames of cached KEGG and mapping files
#' @importFrom BiocFileCache BiocFileCache bfcinfo
#' @export
return_all_cached <- function() {
  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)
  cache_info <- list(
    kegg = BiocFileCache::bfcinfo(bfc_kegg),
    mappings = BiocFileCache::bfcinfo(bfc_map)
  )
  return(cache_info)
}

#' Reset KEGG and mapping caches by deleting all cached files
#' @return None
#' @importFrom BiocFileCache BiocFileCache bfcinfo bfcremove
#' @export
reset_cache <- function() {
  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)
  if (nrow(BiocFileCache::bfcinfo(bfc_kegg)) == 0 && nrow(BiocFileCache::bfcinfo(bfc_map)) == 0) {
    message("No cached files found. Nothing to delete.")
    return()
  }
  message("Deleting all cached KEGG files...")
  message("total kegg pathways files to delete: ", nrow(BiocFileCache::bfcinfo(bfc_kegg)))
  message("total other files to delete: ", nrow(BiocFileCache::bfcinfo(bfc_map)))
  askYesNo("Are you sure you want to delete all cached files?") -> answer
  if (!answer) {
    message("Cache reset aborted.")
    return()
  }
  BiocFileCache::bfcremove(bfc_kegg, BiocFileCache::bfcinfo(bfc_kegg)$rid)
  BiocFileCache::bfcremove(bfc_map, BiocFileCache::bfcinfo(bfc_map)$rid)
}
