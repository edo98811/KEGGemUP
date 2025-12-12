#' Convert KEGG pathway ID to readable pathway name
#'
#' @param id KEGG pathway ID (e.g., 'hsa04110')
#' @return A character string with the readable pathway name
#' @importFrom KEGGREST keggGet
#' @noRd
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
#' @param pathway_id KEGG pathway ID (e.g., 'hsa04110').
#' @param bfc BiocFileCache object for caching KEGG KGML files.
#' @param directory Optional directory to save the KGML file if not using cache.
#' @return Path to the cached KGML file.
#'
#' @importFrom httr2 request req_perform resp_status resp_body_xml resp_is_error req_retry
#' @importFrom BiocFileCache bfcquery bfcpath bfcadd
#' @importFrom xml2 write_xml
#'
#' @examples
#' data_dir <- tempdir()
#' kgml_path <- download_kgml("hsa04110", file_path = data_dir)
#' @export
download_kgml <- function(pathway_id, bfc = NULL, directory = NULL) {

  # check input validity
  if (!is.null(bfc) && !is.null(directory)) {
    stop("Provide either 'bfc' OR 'directory', not both.")
  } else if (!is.null(bfc)) {
    # Check that bfc is a BiocFileCache object
    if (!inherits(bfc, "BiocFileCache")) {
      stop("'bfc' must be a valid BiocFileCache object.")
    }
    mode <- "cache"
  } else if (!is.null(directory)) {
    # Check that directory is a single string
    if (!is.character(directory) || length(directory) != 1) {
      stop("'directory' must be a single string specifying a valid path.")
    }
    # Optionally, create the directory if it does not exist
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
      message("Created directory: ", directory)
    }
    mode <- "dir"
  } else {
    stop("Either 'directory' or 'bfc' must be provided.")
  }

  # Cache key / name
  rname <- paste0(pathway_id, ".xml")

  if (mode == "cache") {
    # Check cache
    qr <- bfcquery(bfc, rname, field = "rname")

    # Check if file exists
    if (nrow(qr) > 0) {
      cached_path <- bfcpath(bfc, qr$rid[1])

      if (file.exists(cached_path)) {
        # if file exists return path
        message("Using cached KEGG KGML for ", pathway_id)
        return(cached_path)
      } else {
        # file missing, re-download
        message("Cache entry found but file missing. Re-downloading.")
      }
    }

    # Temporary file to store KGML
    file_name <- tempfile(fileext = ".xml")
  } else {
    file_name <- path.expand(file.path(directory, rname)) # https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/path.expand
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
      "Failed to download KGML from URL: ", url, " (HTTP status ", resp_status(resp),
      ")"
    )
    return(NULL)
  }

  # Get content as raw vector and check
  kgml_xml <- resp_body_xml(resp)

  write_xml(kgml_xml, file_name)

  if (mode == "cache") {
    # add to BiocFileCache
    res <- bfcadd(bfc, rname = rname, fpath = file_name, action = "copy")
    rid <- names(res)
    message("Downloaded & cached: ", pathway_id)
    return(bfcpath(bfc, rid))
  } else {
    # else return path
    return(file_name)
  }
}

#' Get KEGG compounds with caching.
#' @param bfc A BiocFileCache object for caching.
#' @return A named character vector of KEGG compounds.
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew
#' @importFrom KEGGREST keggList
#' @noRd
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
#' @noRd
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
