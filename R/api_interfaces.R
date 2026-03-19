#' Convert KEGG pathway ID to readable pathway name
#'
#' @param id KEGG pathway ID (e.g., 'hsa04110')
#' @param verbose Logical, if TRUE, print additional messages.
#'
#' @return A character string with the readable pathway name
#'
#' @importFrom KEGGREST keggGet
#'
#' @noRd
get_pathway_name <- function(id, verbose = FALSE) {
  tryCatch(
    {
      res <- KEGGREST::keggGet(id)
      if (verbose) message("Retrieved pathway name for ID: ", id)
      res[[1]]$NAME
    },
    error = function(e) {
      warning(paste("Could not retrieve pathway name for ID:", id))
      return("")
    }
  )
}

#' Download and cache KEGG KGML files.
#'
#' @param pathway_id KEGG pathway ID (e.g., 'hsa04110').
#' @param bfc BiocFileCache object for caching KEGG KGML files.
#' @param directory Optional directory to save the KGML file if not using cache.
#' @param verbose Logical, if TRUE, print additional messages.
#'
#' @return Path to the cached KGML file.
#'
#' @importFrom httr2 request req_perform resp_status resp_body_xml resp_is_error req_retry
#' @importFrom BiocFileCache bfcquery bfcpath bfcadd
#' @importFrom xml2 write_xml
#'
#' @export
#'
#' @examples
#' data_dir <- tempdir()
#' kgml_path <- download_kgml("hsa04110", directory = data_dir, verbose = TRUE)
download_kgml <- function(pathway_id, bfc = NULL, directory = NULL, verbose = FALSE) {
  # check input validity
  if (!is.null(bfc) && !is.null(directory)) {
    stop("Provide either 'bfc' OR 'directory', not both.")
  } else if (!is.null(bfc)) {
    if (!inherits(bfc, "BiocFileCache")) {
      stop("'bfc' must be a valid BiocFileCache object.")
    }
    mode <- "cache"
  } else if (!is.null(directory)) {
    if (!is.character(directory) || length(directory) != 1) {
      stop("'directory' must be a single string specifying a valid path.")
    }
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
      if (verbose) message("Created directory: ", directory)
    }
    mode <- "dir"
  } else {
    directory <- getwd()
    if (verbose) message("No 'bfc' or 'directory' provided. Using current working directory: ", directory)
    mode <- "dir"
  }

  if (!is_valid_pathway(pathway_id)) {
    stop("Invalid KEGG pathway ID format.")
  }

  url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
  if (verbose) message("Downloading KGML from: ", url)

  if (mode == "cache") {
    path <- BiocFileCache::bfcrpath(bfc, url, ext = ".xml")
    if (verbose) message("Cached: ", pathway_id)
    return(path)
  } else {
    file_name <-
      if (grepl("\\.[^/\\\\]+$", directory)) {
        directory
      } else {
        file.path(path.expand(directory), paste0(pathway_id, ".xml"))
      }
    resp <- request(url) |>
      req_retry(max_tries = 3) |>
      req_perform()

    if (resp_is_error(resp)) {
      warning(
        "Failed to download KGML from URL: ", url, " (HTTP status ", resp_status(resp),
        ")"
      )
      return(NULL)
    }

    kgml_xml <- resp_body_xml(resp)
    write_xml(kgml_xml, file_name)
    if (verbose) message("Downloaded & saved in: ", file_name)

    return(file_name)
  }
}

#' Get KEGG db with caching.
#'
#' @param db_name Name of the KEGG database to retrieve (default: "compound").
#' @param bfc BiocFileCache object for caching KEGG database files.
#' @param directory Optional directory to save the KEGG database file if not using cache.
#' @param verbose Logical, if TRUE, print additional messages.
#' @return A data frame with KEGG IDs and names.
#' @details The valid KEGG database names are:
#' kegg | pathway | brite | module | ko | genes | <org> | vg | vp | ag |
#' genome | ligand | compound | glycan | reaction | rclass | enzyme |
#' network | variant | disease | drug | dgroup
#' @importFrom KEGGREST keggList
#' @importFrom utils read.table write.table
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcpath bfcnew bfcadd bfcrpath
#' @importFrom httr2 request req_perform resp_status resp_body_string resp_is_error req_retry
#' @examples
#' # Saving in directory
#' data_dir <- tempdir()
#' kegg_compounds <- get_kegg_db("compound", directory = data_dir, verbose = TRUE)
#' # Just returning without saving
#' kegg_genes <- get_kegg_db("compound", verbose = TRUE)
#' @export
get_kegg_db <- function(
    db_name = "compound",
    directory = NULL,
    bfc = NULL,
    verbose = FALSE) {
  if (verbose) message("Retrieving KEGG database: ", db_name)
  if (!is.null(bfc) && !is.null(directory)) {
    stop("Provide either 'bfc' OR 'directory', not both.")
  } else if (!is.null(bfc)) {
    if (!inherits(bfc, "BiocFileCache")) {
      stop("'bfc' must be a valid BiocFileCache object.")
    }
    mode <- "cache"
  } else if (!is.null(directory)) {
    if (!is.character(directory) || length(directory) != 1) {
      stop("'directory' must be a single string specifying a valid path.")
    }
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
      if (verbose) message("Created directory: ", directory)
    }
    mode <- "dir"
  } else {
    if (verbose) message("No 'bfc' or 'directory' provided. Not saving KEGG database only downloading and returning.")
    mode <- "none"
  }

  url <- paste0("https://rest.kegg.jp/list/", db_name)

  if (mode == "cache") {
    path <- BiocFileCache::bfcrpath(bfc, url, ext = ".tsv")
    if (verbose) message("Cached KEGG database: ", db_name)
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
    if (verbose) message("Downloaded & saved KEGG database in: ", file_name)
  }

  return(kegg_db)
}

#' Download all pathways
#'
#' Download all KEGG pathways for a given organism
#'
#' @details
#' This function can take a while to run completely, but is an efficient way
#' to retrieve all the kgml files encoding for the KEGG pathways.
#' Some failures could be triggered by a rate limitation imposed by the KEGG
#' website. Since this caches the files locally, it is safe to rerun to complete
#' the operation without extra unneeded requests to the API
#'
#' @param org KEGG organism code (e.g., 'hsa' for human).
#' @param verbose Logical, if TRUE, print additional messages.
#'
#' @return The BiocFileCache object is returned invisibly
#'
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom utils askYesNo
#'
#' @export
#'
#' @examples
#' # Download all pathways for human
#' # download_all_pathways("hsa", verbose = TRUE)
download_all_pathways <- function(org, verbose = FALSE) {
  path <- tools::R_user_dir("BiocFileCache", which = "cache")
  bfc_kegg <- BiocFileCache(cache = file.path(path, "kegg_maps"), ask = FALSE)
  bfc_map <- BiocFileCache(cache = file.path(path, "mappings"), ask = FALSE)

  all_pathways_df <- get_kegg_db(bfc = bfc_map,
                                 db_name = paste0("pathway/", org),
                                 verbose = verbose)

  answer <- askYesNo(
    msg = paste0("Download all ", nrow(all_pathways_df),
                 " KEGG pathways for organism '", org, "'? This may take a while."))

  if (!answer) {
    if (verbose) message("Aborting download of all pathways.")
    return(NULL)
  }

  tot_pathways <- nrow(all_pathways_df)

  for (i in seq_len(tot_pathways)) {
    pathway_id <- all_pathways_df$kegg_id[i]
    pathway_desc <- all_pathways_df$description[i]
    if (verbose) message(i, "/", tot_pathways, " - ", pathway_id, "|", pathway_desc)
    download_kgml(pathway_id, bfc = bfc_kegg, verbose = verbose)
  }

  message("Done retrieving all pathways for ", org, "!")

  # return invisibly the BFC object, enabling further processing
  invisible(bfc_kegg)
}
