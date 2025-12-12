#' Check the structure and validity of each entry in de_results
#'
#' @param de_entry An entry from the de_results list
#' @param name The name of the entry in the de_results list (for error messages)
#'
#' @return Invisible TRUE if all checks pass, otherwise stops with an error
is_valid_de_entry <- function(de_entry, name) {
  # --- Structure checks ---
  if (!is.list(de_entry) || !all(c("de_table", "value_column", "feature_column") %in%
    names(de_entry))) {
    warning(
      "Each entry in de_results must be a list with elements:
      de_table, value_column, feature_column (problem in '",
      name, "')"
    )
    return(FALSE)
  }
  if (!inherits(de_entry$de_table, "data.frame")) {
    warning("de_results[['", name, "']]$de_table must be a data frame")
    return(FALSE)
  }
  # Check that value_column is character and present in de_table
  if (!is.character(de_entry$value_column) || !(de_entry$value_column %in% colnames(de_entry$de_table))) {
    warning(
      "de_results[['", name, "']]$value_column must be a column name in de_results[['",
      name, "']]$de_table"
    )
    return(FALSE)
  }
  # Check feature_column is character and present in de_table or is
  # 'rownames'
  if (!is.character(de_entry$feature_column) || !(de_entry$feature_column %in%
    c(colnames(de_entry$de_table), "rownames"))) {
    warning(
      "de_results[['", name, "']]$feature_column must be a column name in de_results[['",
      name, "']]$de_table or 'rownames'"
    )
    return(FALSE)
  }

  # --- Row check ---

  # if feature_column is not in the colnames, if it is null, if it is not
  # 'rownames' then error
  if (!(de_entry$feature_column %in% colnames(de_entry$de_table)) && de_entry$feature_column !=
    "rownames") {
    warning(paste(
      "Column", de_entry$feature_column, "not found in de_results: ",
      name
    ))
    return(FALSE)
  }

  # --- Check that table is correct ---
  if (is.null(de_entry$de_table) || !(de_entry$value_column %in% colnames(de_entry$de_table))) {
    warning("Invalid de_table or value_column in de_results: ", name)
    return(FALSE)
  }

  return(TRUE)
}

#' Validate KGML file structure
#' @param file_path Path to the KGML file
#' @return TRUE if valid KGML, FALSE otherwise
is_valid_kgml <- function(file_path) {
  # Try to read the XML file
  xml_content <- tryCatch(
    {
      xml2::read_xml(file_path)
    },
    error = function(e) {
      warning(paste("Failed to read KGML file:", file_path))
      return(NULL)
    }
  )

  if (is.null(xml_content)) {
    return(FALSE)
  }

  # Check for the presence of the root 'pathway' node
  root_node <- xml2::xml_find_first(xml_content, "/pathway")
  if (is.na(root_node)) {
    warning(paste(
      "KGML file does not contain a valid 'pathway' root node:",
      file_path
    ))
    return(FALSE)
  }

  return(TRUE)
}

#' Validate KEGG pathway ID format
#' @param pathway_id KEGG pathway ID (e.g., 'hsa04110')
#' @return TRUE if valid format, FALSE otherwise
is_valid_pathway <- function(pathway_id) {
  # Check if pathway_id matches KEGG pathway formats: 'hsa04110' or '04110'
  if (!is.character(pathway_id) || length(pathway_id) != 1) {
    return(FALSE)
  }
  grepl("^[a-z]{2,3}\\d{5}$", pathway_id)
}
