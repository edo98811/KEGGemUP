#' Convert KEGG IDs with prefixes to a single string of IDs without prefixes
#' @param kegg_ids Character string of KEGG IDs with prefixes
#' (e.g., 'cpd:C00022 mmu:1234 ko:K00001 ko:K00002')
#' @return Character string of KEGG IDs without prefixes, separated by ';'
#' @noRd
remove_kegg_prefix_str <- function(kegg_ids) {
  sapply(kegg_ids, function(id) {
    if (is.na(id)) {
      return(NA_character_)
    } # preserve NA
    elements <- strsplit(id, " ")[[1]] # split by space
    elements <- sub("^[a-z]+:", "", elements) # remove prefix
    paste(elements, collapse = ";") # collapse back to single string
  }, USE.NAMES = FALSE)
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
#' @details This function retrieves information about
#' all cached KEGG pathway files
#' and mapping files stored using BiocFileCache.
#' @examples
#' cache_info <- return_all_cached()
#' print(cache_info$kegg) # View cached KEGG pathway files
#' print(cache_info$mappings) # View cached mapping files
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
#' Create a mapping data frame from vertices_df for subgraph extraction
#' @param vertices_df Data frame of graph vertices with columns 'name' and 'ids_for_mapping'
#' @return A data frame with columns 'name' and 'matched_id'. 
#' The 'matched_id' column contains individual KEGG IDs extracted from 'ids_for_mapping', which may contain multiple IDs separated by ';'.
#' @details Example of output:
#'   name matched_id
#' 1  NodeA     K00001
#' 2  NodeA     K00002 
#' 3  NodeB     K00003
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

#' Reset KEGG and mapping caches by deleting all cached files
#' @return None
#' @importFrom BiocFileCache BiocFileCache bfcinfo bfcremove
#' @importFrom utils askYesNo
#' @details This function deletes all cached KEGG pathway files
#' and mapping files stored using BiocFileCache.
#' It prompts the user for confirmation before proceeding with the deletion.
#' @examples
#' # reset_cache()
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
  message("total kegg pathways files to delete: ", paste(BiocFileCache::bfcinfo(bfc_kegg)$rname, sep = ", "))
  message("total other files to delete: ", paste(BiocFileCache::bfcinfo(bfc_map)$rname, sep = ", "))
  askYesNo("Are you sure you want to delete all cached files?") -> answer
  if (!answer) {
    message("Cache reset aborted.")
    return(invisible(NULL))
  }
  BiocFileCache::bfcremove(bfc_kegg, BiocFileCache::bfcinfo(bfc_kegg)$rid)
  BiocFileCache::bfcremove(bfc_map, BiocFileCache::bfcinfo(bfc_map)$rid)
}
