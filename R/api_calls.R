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

# #' Download KEGG KGML file for a given pathway ID, with caching.
# #'
# #' @param url URL to download the KGML file from.
# #' @return Content of the respose as text, or NULL if download failed.
# #' @importFrom httr GET http_error content status_code
# get_kgml <- function(url) {
#   # Download
#   resp <- httr::GET(url)

#   # Check if the request was successful (status 200)
#   if (httr::http_error(resp)) {
#     warning(
#       "Failed to download KGML from URL: ", url,
#       " (HTTP status ", httr::status_code(resp), ")"
#     )
#     return(NULL)
#   }

#   content <- httr::content(resp, as = "raw", encoding = "")
#   return(content)
# }


#' Download and cache KEGG KGML files.
#' @param pathway_id KEGG pathway ID (e.g., "hsa04110").
#' @param bfc BiocFileCache object for caching KEGG KGML files.
#' @return Path to the cached KGML file.
#'
#' @importFrom xml2 read_xml
#' @importFrom httr GET http_error content status_code
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew bfcadd
get_and_cache_kgml <- function(pathway_id, bfc) {
  # cache key / name
  rname <- paste0(pathway_id, ".xml")

  # Check cache
  qr <- bfcquery(bfc, rname, field = "rname")

  if (nrow(qr) > 0) {
    message("Using cached KEGG KGML for ", pathway_id)
    return(bfcpath(bfc, qr$rid[1]))
  }

  # If not cached, download
  message("Downloading KGML for ", pathway_id, "...")

  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")

  # # Function to get KGML content
  # kgml_content <- get_kgml(url)

  # Make request
  resp <- httr::GET(url)

  # Check if the request was successful (status 200)
  if (httr::http_error(resp)) {
    warning(
      "Failed to download KGML from URL: ", url,
      " (HTTP status ", httr::status_code(resp), ")"
    )
    return(NULL)
  }

  # Get content as raw vector and check
  kgml_raw <- httr::content(resp, as = "raw", encoding = "")
  if (is.null(kgml_raw)) {
    warning("Failed to download KGML for ", pathway_id)
    return(NULL)
  }

  # # create cache entry
  # dest <- bfcnew(bfc, rname = rname, ext = ".xml")

  # # write content
  # writeBin(kgml_content, dest)

  # ---- WINDOWS SAFE WAY ----
  # Write to a temp file IN BINARY MODE (Windows-safe)
  tmp <- tempfile(fileext = ".xml")
  con <- file(tmp, "wb")
  writeBin(kgml_raw, con)
  close(con)

  # Import into BiocFileCache
  res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "move")
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