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
#' @importFrom httr GET http_error content status_code
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew bfcadd
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

  txt <- rawToChar(kgml_raw)
  tmp <- tempfile(fileext = ".xml")

  writeLines(txt, con = tmp)

  # Add to BiocFileCache
  res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "copy")
  rid <- names(res)

  message("Downloaded & cached: ", pathway_id)
  return(bfcpath(bfc, rid))
}
# get_and_cache_kgml <- function(pathway_id, bfc) {
#   # Cache key / name
#   rname <- paste0(pathway_id, ".xml")

#   # Check cache
#   qr <- bfcquery(bfc, rname, field = "rname")
#   if (nrow(qr) > 0) {
#     message("Using cached KEGG KGML for ", pathway_id)
#     return(bfcpath(bfc, qr$rid[1]))
#   }

#   # Download KGML
#   message("Downloading KGML for ", pathway_id, " ...")
#   url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
#   tmp <- tempfile(fileext = ".xml")
#   download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)

#   # Add to BiocFileCache
#   res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "copy")
#   rid <- names(res)

#   message("Downloaded & cached: ", pathway_id)
#   return(bfcpath(bfc, rid))
# }

# get_and_cache_kgml <- function(pathway_id, bfc) {
#   # cache key / name
#   rname <- paste0(pathway_id, ".xml")

#   # Check cache
#   qr <- bfcquery(bfc, rname, field = "rname")

#   if (nrow(qr) > 0) {
#     message("Using cached KEGG KGML for ", pathway_id)
#     return(bfcpath(bfc, qr$rid[1]))
#   }

#   # If not cached, download
#   message("Downloading KGML for ", pathway_id, "...")

#   url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")

#   # # Function to get KGML content
#   # kgml_content <- get_kgml(url)

#   # Make request
#   resp <- httr::GET(url)

#   # Check if the request was successful (status 200)
#   if (httr::http_error(resp)) {
#     warning(
#       "Failed to download KGML from URL: ", url,
#       " (HTTP status ", httr::status_code(resp), ")"
#     )
#     return(NULL)
#   }

#   # Get content as raw vector and check
#   kgml_raw <- httr::content(resp, as = "raw", encoding = "")
#   if (is.null(kgml_raw)) {
#     warning("Failed to download KGML for ", pathway_id)
#     return(NULL)
#   }

  # # Write raw content to temporary file
  # con <- rawConnection(kgml_raw)
  # xml_lines <- readLines(con, encoding = "UTF-8")
  # close(con)
  # kgml_text <- rawToChar(kgml_raw)
  # kgml_clean <- sub("\\n$", "", kgml_text) # remove trailing newline
  # kgml_clean <- sub('\\0', '', kgml_clean)
  # kgml_clean_raw <- charToRaw(kgml_clean)

  # xml_lines <- head(xml_lines, -1) # remove empty last line (alaways present)

  # Write the XML lines directly
#   txt <- rawToChar(kgml_raw)
#   end_tag <- "</pathway>"
#   pos <- regexpr(end_tag, txt, fixed = TRUE)
#   if (pos[1] == -1) stop("No </pathway> tag found")
#   truncated <- substr(txt, 1, pos[1] + attr(pos, "match.length") - 1)
#   # truncated <- sub("[\r\n\x00\\s]+$", "", truncated, perl = TRUE)
#   kgml_clean <- charToRaw(truncated)

#   kgml_clean <- kgml_clean[kgml_clean != as.raw(0)]
#   tmp <- tempfile(fileext = ".xml")
#   con <- file(tmp, "wb")
#   writeBin(kgml_clean, con)
#   close(con)

#   # Add to cache
#   res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "move")
#   rid <- names(res)

#   message("Downloaded & cached: ", pathway_id)
#   return(bfcpath(bfc, rid))
# }
# return(dest)

# dest <- bfcnew(bfc, rname = rname, ext = ".xml")
# # create cache entry
# dest <- bfcnew(bfc, rname = rname, ext = ".xml")

# # write content
# writeBin(kgml_content, dest)

# ---- WINDOWS SAFE WAY ----
# Write to a temp file IN BINARY MODE (Windows-safe)
# kgml_text <- rawToChar(kgml_raw)
# Encoding(kgml_text) <- "UTF-8" # read or Set the Declared Encodings for a Character Vector
# kgml_text <- sub("(?s)(</pathway>).*", "\\1", kgml_text, perl = TRUE) # remove anything after the last closing tag
# kgml_clean <- charToRaw(kgml_text) # write back to binary
# # tmp <- tempfile(fileext = ".xml")
# # con <- file(tmp, "wb")
# # writeBin(kgml_raw, con)
# # close(con)
# dest <- bfcnew(bfc, rname = rname, ext = ".xml")
# writeBin(kgml_clean, dest)
# # # Import into BiocFileCache
# # res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "move")
# # rid <- names(res)

# message("Downloaded & cached: ", pathway_id)
# # return(bfcpath(bfc, rid))
# return(dest)
# Find the last </pathway> tag in raw bytes
# close_tag <- charToRaw("</pathway>")
# pos <- max(gregexpr(close_tag, kgml_raw, fixed = TRUE)[[1]])
# if (pos > 0) {
#   kgml_clean <- kgml_raw[1:(pos + length(close_tag) - 1)]
# } else {
#   kgml_clean <- kgml_raw
# }

# # Write to temporary file
# tmp <- tempfile(fileext = ".xml")
# conn <- file(tmp, "wb")
# writeBin(kgml_clean, conn)
# close(conn)

# Add to cache safely
# res <- bfcadd(bfc, rname = rname, fpath = tmp, action = "move")
# rid <- names(res)

# message("Downloaded & cached: ", pathway_id)
# return(bfcpath(bfc, rid))
# return(dest)

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
