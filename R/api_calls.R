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
#' kgml_path <- download_kgml("hsa04110", directory = data_dir)
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

  # Validate pathway ID format
  if (!is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
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

#' Get KEGG db with caching.
#'
#' @param bfc A BiocFileCache object for caching.
#' @param db_name KEGG database name (e.g., 'compound', 'glycan').
#' @return A data frame with KEGG IDs and names.
#' @details The valid KEGG database names are:  
#' kegg | pathway | brite | module | ko | genes | <org> | vg | vp | ag |
#' genome | ligand | compound | glycan | reaction | rclass | enzyme |
#' network | variant | disease | drug | dgroup
#' @importFrom KEGGREST keggList
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew bfcadd bfcrpath
#' @noRd
get_kegg_db <- function(bfc, db_name = "compound") {

  url <- paste0("https://rest.kegg.jp/list/", db_name)

  cache_name <- paste0(db_name, ".rds")
  
  path <- BiocFileCache::bfcrpath(bfc, url)
  kegg_db <- read.table(path, sep = "\t") |> data.frame()    

  # # Check if cache exists
  # qr <- BiocFileCache::bfcquery(bfc, cache_name, field = "rname")
  # url <- paste0("https://rest.kegg.jp/list/", db_name)

  # if (nrow(qr) > 0) {
  #   message("Loading KEGG ", db_name, " from cache...")
  #   kegg_db <- readRDS(BiocFileCache::bfcpath(bfc, qr$rid[1]))
  #   return(kegg_db)
  # }

  # # Otherwise download from KEGG
  # message("Downloading KEGG ", db_name, "...")
  # kegg_db <- KEGGREST::keggList(db_name)

  # temp_file_path <- tempfile(fileext = ".rds")
  # saveRDS(kegg_db, file = temp_file_path)

  # # Save to cache
  # res <- bfcadd(bfc, rname = cache_name, fpath = temp_file_path, action = "copy")

  return(kegg_db)
}