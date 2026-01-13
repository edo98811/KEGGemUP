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
    directory <- getwd()
    message("No 'bfc' or 'directory' provided. Using current working directory: ", directory)
    mode <- "dir"
  }

  # Validate pathway ID format
  if (!is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
  }

  # Cache key / name
  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
  message("Downloading KGML from: ", url)

  if (mode == "cache") {
    path <- BiocFileCache::bfcrpath(bfc, url, ext = ".xml")

    message("Downloaded & cached: ", pathway_id)
    return(path)
  } else {
    file_name <-
      if (grepl("\\.[^/\\\\]+$", directory)) {
        directory
      } else {
        file.path(path.expand(directory), paste0(pathway_id, ".xml")) # https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/path.expand
      }
    resp <- request(url) |>
      req_retry(max_tries = 3) |>
      req_perform()

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
    message("Downloaded & saved in: ", file_name)

    return(file_name)
  }
}

#' Get KEGG db with caching.
#'
#' @param db_name Name of the KEGG database to retrieve (default: "compound").
#' @param bfc BiocFileCache object for caching KEGG database files.
#' @param directory Optional directory to save the KEGG database file if not using cache.
#' @return A data frame with KEGG IDs and names.
#' @details The valid KEGG database names are:
#' kegg | pathway | brite | module | ko | genes | <org> | vg | vp | ag |
#' genome | ligand | compound | glycan | reaction | rclass | enzyme |
#' network | variant | disease | drug | dgroup
#' @details If neither 'bfc' nor 'directory' is provided, the KEGG database
#' will be downloaded but not saved. It will be returned as a data frame.
#' @importFrom KEGGREST keggList
#' @importFrom utils read.table write.table
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew bfcadd bfcrpath
#' @importFrom httr2 request req_perform resp_status resp_body_string resp_is_error req_retry
#' @export
get_kegg_db <- function(db_name = "compound", directory = NULL, bfc = NULL) {
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
    # Create the directory if it does not exist
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
      message("Created directory: ", directory)
    }
    mode <- "dir"
  } else {
    message("No 'bfc' or 'directory' provided. Not saving KEGG database only downloading and returning.")
    mode <- "none"
  }

  url <- paste0("https://rest.kegg.jp/list/", db_name)

  if (mode == "cache") {
    path <- BiocFileCache::bfcrpath(bfc, url, ext = ".tsv")
    message("Downloaded & cached KEGG database: ", db_name)
    con <- path
  } else {
    resp <- request(url) |>
      req_retry(max_tries = 3) |>
      req_perform()

    if (resp_is_error(resp)) {
      warning(
        "Failed to download KEGG DB: ", db_name,
        " (HTTP status ", resp_status(resp), ")"
      )
      return(NULL)
    }
    con <- textConnection(httr2::resp_body_string(resp))
    on.exit(close(con), add = TRUE)
  }

  kegg_db <- read.table(
    con,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    col.names = c("kegg_id", "description")
  ) |> as.data.frame()

  if (mode == "dir") {
    file_name <-
      if (grepl("\\.[^/\\\\]+$", directory)) {
        directory
      } else {
        file.path(path.expand(directory), paste0("kegg_", db_name, ".tsv"))
      }
    write.table(
      kegg_db,
      file = file_name,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    message("Downloaded & saved KEGG database in: ", file_name)
  }

  return(kegg_db)
}

#' Download all KEGG pathways for a given organism.
#' @param org KEGG organism code (e.g., 'hsa' for human).
#' @return None
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom utils askYesNo
#' @export
download_all_pathways <- function(org) {
  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)

  all_pathways <- get_kegg_db(bfc_map, paste0("pathway/", org))

  askYesNo("Download all ", nrow(all_pathways), "KEGG pathways for organism '", org, "'? This may take a while.") -> answer
  if (!answer) {
    message("Aborting download of all pathways.")
    return(NULL)
  }

  for (pathway_id in all_pathways[, 1]) {
    download_kgml(pathway_id, bfc = bfc_kegg)
  }
}
