#' Convert KEGG pathway ID to readable pathway name
#'
#' @param id KEGG pathway ID (e.g., "hsa04110")
#' @return A character string with the readable pathway name
get_pathway_name <- function(id) {
  tryCatch(
    {
      res <- KEGGREST::keggGet(id)
      res[[1]]$NAME
    },
    error = function(e) {
      warning(paste("Could not retrieve pathway name for ID:", id))
      return("")
    }
  )
}

#' Download and cache KEGG KGML files.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param bfc BiocFileCache object for caching KEGG KGML files.
#' @return Path to the cached KGML file.
#'
#' @importFrom httr2 request req_perform resp_status resp_body_xml resp_is_error req_retry
#' @importFrom BiocFileCache bfcquery bfcpath bfcadd
#' @importFrom xml2 write_xml
get_and_cache_kgml <- function(pathway_id, bfc) {
  # Cache key / name
  rname <- paste0(pathway_id, ".xml")

  # Check cache
  qr <- bfcquery(bfc, rname, field = "rname")
  if (nrow(qr) > 0) {
    message("Using cached KEGG KGML for ", pathway_id)
    return(bfcpath(bfc, qr$rid[1]))
  }

  # Download KGML
  message("Downloading KGML for ", pathway_id, " ...")
  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
  resp <- request(url) |>
    req_retry(max_tries = 3) |>
    req_perform(error_call = FALSE)

  # Check success
  if (resp_is_error(resp)) {
    warning(
      "Failed to download KGML from URL: ", url,
      " (HTTP status ", resp_status(resp), ")"
    )
    return(NULL)
  }

  # Get content as raw vector and check
  kgml_xml <- resp_body_xml(resp)

  # Optionally: write to cache file
  tmp <- tempfile(fileext = ".xml")
  write_xml(kgml_xml, tmp) # preserves formatting

  # Add to BiocFileCache
  res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "copy")
  rid <- names(res)

  message("Downloaded & cached: ", pathway_id)
  return(bfcpath(bfc, rid))
}

#' Get KEGG compounds with caching.
#' @param bfc A BiocFileCache object for caching.
#' @return A named character vector of KEGG compounds.
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew
#' @importFrom KEGGREST keggList
get_kegg_compounds <- function(bfc) {
  cache_name <- "kegg_compounds.rds"

  # Check if cache exists
  qr <- bfcquery(bfc, cache_name, field = "rname")

  if (nrow(qr) > 0) {
    message("Loading KEGG compounds from cache...")
    compounds <- readRDS(bfcpath(bfc, qr$rid[1]))
    return(compounds)
  }

  # Otherwise download from KEGG
  message("Downloading KEGG compounds...")
  kegg_compounds <- KEGGREST::keggList("compound")

  # Save to cache
  path <- bfcnew(bfc, rname = cache_name, ext = ".rds")
  saveRDS(kegg_compounds, file = path)

  return(kegg_compounds)
}

#' Get KEGG glycans with caching.
#'
#' @param bfc A BiocFileCache object for caching.
#' @return A data frame with KEGG glycan IDs and names.
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew
#' @importFrom KEGGREST keggList
get_kegg_glycans <- function(bfc) {
  cache_name <- "kegg_glycans.rds"

  # Check if cache exists
  qr <- BiocFileCache::bfcquery(bfc, cache_name, field = "rname")

  if (nrow(qr) > 0) {
    message("Loading KEGG glycans from cache...")
    glycans <- readRDS(BiocFileCache::bfcpath(bfc, qr$rid[1]))
    return(glycans)
  }

  # Otherwise download from KEGG
  message("Downloading KEGG glycans...")
  kegg_glycans <- KEGGREST::keggList("glycan")

  # Save to cache
  path <- BiocFileCache::bfcnew(bfc, rname = cache_name, ext = ".rds")
  saveRDS(kegg_glycans, file = path)

  return(kegg_glycans)
}
