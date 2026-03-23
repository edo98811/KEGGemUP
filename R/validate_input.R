#' Check the structure and validity of each entry in de_results
#'
#' @param de_entry An entry from the de_results list
#' @param name The name of the entry in the de_results list (for error messages)
#'
#' @return TRUE if all checks pass, otherwise stops with an error
#' @noRd
is_valid_de_entry <- function(de_entry, name) {
  #  Structure checks
  if (!is.list(de_entry) || !all(c("de_table", "value_column", "feature_column") %in%
    names(de_entry))) {
    warning(
      "Each entry in de_results must be a list with elements:
      de_table, value_column, feature_column (problem in '",
      name, "')"
    )
    return(FALSE)
  }

  if (!is_valid_dataframe(de_entry$de_table, name)) {
    return(FALSE)
  }

  # Check feature_column is character and present in de_table or is
  # 'rownames'
  if (!is.character(de_entry$feature_column) || !(de_entry$feature_column %in%
    c(colnames(de_entry$de_table), "rownames"))) {
    warning(
      "de_results list element must have a valid feature_column (problem in '", name, "')"
    )
    return(FALSE)
  }

  #  Row chec
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

  #  Check that table is correct
  if (is.null(de_entry$de_table) || !(de_entry$value_column %in% colnames(de_entry$de_table))) {
    warning("Invalid de_table or value_column in de_results: ", name)
    return(FALSE)
  }

  return(TRUE)
}

#' Validate KGML file structure
#' @param file_path Path to the KGML file
#' @return TRUE if valid KGML, FALSE otherwise
#' @importFrom xml2 read_xml xml_find_first
#' @noRd
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
#' @noRd
is_valid_pathway <- function(pathway_id) {
  # Check if pathway_id matches KEGG pathway formats: 'hsa04110' or '04110'
  if (!is.character(pathway_id) || length(pathway_id) != 1) {
    return(FALSE)
  }
  grepl("^[a-z]{2,3}\\d{5}$", pathway_id)
}


is_valid_dataframe <- function(obj, name = "object") {
  # Non-null
  if (is.null(obj)) {
    warning(name, " is NULL")
    return(FALSE)
  }

  # Single object (length 1 for atomic, for list/data.frame it’s rows > 0)
  if (is.atomic(obj) && length(obj) != 1) {
    warning(name, " is not a single object, make sure you are not providing MArrayLM but the uptput of topTable, if using limma")
    return(FALSE)
  }

  # Must inherit from data.frame
  if (!inherits(obj, "data.frame")) {
    warning(name, " must be a data.frame, but is of class: ", paste(class(obj), collapse = "/"))
    return(FALSE)
  }

  # Non-empty (has rows and columns)
  if (nrow(obj) == 0 || ncol(obj) == 0) {
    warning(name, " is empty (0 rows or 0 columns)")
    return(FALSE)
  }

  # Column names exist
  if (is.null(colnames(obj)) || any(colnames(obj) == "")) {
    warning(name, " has missing or empty column names")
    return(FALSE)
  }

  # Passed all checks
  return(TRUE)
}

# standardize_de_results <- function(
#   de_results,
#   value_column = NULL,
#   feature_column = NULL
# ) {
#   # If it's a data.frame, value_column and feature_column must be provided as character strings,
#   # otherwise if it's a list, value_column and feature_column should not be provided
#   # invalid_input <- !(
#   #   xor(
#   #     is.data.frame(de_results) && is.character(value_column) && is.character(feature_column)) &&
#   #     length(value_column) == 1 && length(feature_column) == 1 && !is.null(names(de_results)),
#   #     is.list(de_results) && !is.data.frame(de_results) && is.null(value_column) &&
#   #     is.null(feature_column)
#   # )

#   df_case <-
#     is.data.frame(de_results) &&
#       is.character(value_column) &&
#       is.character(feature_column) &&
#       length(value_column) == 1 &&
#       length(feature_column) == 1

#   list_case <-
#     is.list(de_results) &&
#       !is.data.frame(de_results) &&
#       is.null(value_column) &&
#       is.null(feature_column) &&
#       !is.null(names(de_results))

#   valid_input <- xor(df_case, list_case)

#   if (!valid_input) {
#     warning(
#       "If de_results is a data.frame, value_column and feature_column must be provided as character strings.
#       If de_results is a list, value_column and feature_column should not be provided. Ignoring de_results.
#       The list must be named and each entry must be a list with elements: de_table (data.frame), value_column (character),
#       feature_column (character). Make sure you are not providing MArrayLM but the output of topTable, if using limma."
#     )
#     return(NULL)
#   }

#   # If are single values, check that they are valid and convert to list format
#   if ((is.data.frame(de_results) && is.character(value_column) && is.character(feature_column))) {
#     if (all(
#       data_frame_has_columns(de_results, c(feature_column, value_column)),
#       dataframe_cols_are_of_type(de_results, setNames(
#         c("character", "numeric"),
#         c(feature_column, value_column)
#       ))
#     )) {
#       message(
#         "de_results provided as a single data.frame. ",
#         "Using provided value_column: '", value_column, "' and feature_column: '", feature_column, "'."
#       )

#       de_results <- list(
#         de_input = list(
#           de_table = de_results,
#           value_column = value_column,
#           feature_column = feature_column
#         )
#       )
#     } else {
#       warning("Provided data.frame does not have the required columns with correct types. Ignoring de_results.")
#       return(NULL)
#     }
#   }

#   ## If it's a list, keep only valid entries
#   keep <- vapply(
#     names(de_results),
#     function(name) {
#       de_entry <- de_results[[name]]
#       all(
#         list_has_names(de_entry, c("de_table", "value_column", "feature_column")),
#         list_elements_are_types(de_entry, c(de_table = "data.frame", value_column = "character", feature_column = "character")),
#         data_frame_has_columns(de_entry$de_table, c(de_entry$feature_column, de_entry$value_column)),
#         setNames(
#           c("character", "numeric"),
#           c(de_entry$feature_column, de_entry$value_column)
#         ),
#         length(de_entry$value_column) == 1 && length(de_entry$feature_column) == 1
#       )
#     },
#     logical(1)
#   )

#   de_results <- de_results[keep]

#   # If nothing valid remains, return NULL
#   if (length(de_results) == 0) {
#     return(NULL)
#   }

#   de_results
# }

# ' Normalize de_results input into a standard format
#' @param de_results NULL, a data.frame, or a named list of de_results entries
#' @param value_column Column name for the differential expression values (used if de_results is a data.frame)
#' @param feature_column Column name for the feature IDs (used if de_results is a
#' data.frame)
#' @return A named list of validated de_results entries, or NULL if invalid
#' @noRd
standardize_de_results <- function(
  de_results,
  value_column = NULL,
  feature_column = NULL
) {
  if (is.null(de_results)) {
    return(NULL)
  }

  # data.frame into default named list
  if (is.list(de_results) && !is.data.frame(de_results)) {
    if (is.null(names(de_results)) ||
      any(names(de_results) == "")) {
      warning(
        "de_results must be NULL, a valid data.frame, or a named list. ",
        "Ignoring de_results, make sure you are not providing MArrayLM but the uptput of topTable, if using limma"
      )
      return(NULL)
    }
  } else if (is_valid_dataframe(de_results, name = "de_results")) {
    message(
      "de_results provided as a single data.frame. ",
      if (is.null(value_column)) {
        message("Using default value_column: 'log2FoldChange'")
      } else {
        message(paste0("Using provided value_column: '", value_column, "'"))
      },
      if (is.null(feature_column)) {
        message(" and default feature_column: 'KEGG_ids'.")
      } else {
        message(paste0(" and provided feature_column: '", feature_column, "'."))
      }
    )

    de_results <- list(
      de_input = list(
        de_table = de_results,
        value_column = if (is.null(value_column)) "log2FoldChange" else value_column,
        feature_column = if (is.null(feature_column)) "KEGG_ids" else feature_column
      )
    )
  } else {
    warning(
      "de_results must be NULL, a valid data.frame, or a named list. ",
      "Ignoring de_results."
    )
    return(NULL)
  }

  # keep only valid entries
  keep <- vapply(
    names(de_results),
    function(name) {
      is_valid_de_entry(de_results[[name]], name)
    },
    logical(1)
  )

  de_results <- de_results[keep]

  # If nothing valid remains, return NULL
  if (length(de_results) == 0) {
    return(NULL)
  }

  de_results
}

# Check that a data frame has specific columns of specific types
data_frame_has_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    warning("Data frame is missing required columns: ", paste(missing, collapse = ", "))
    return(FALSE)
  }
  TRUE
}

# # Check that a column has no NA or empty strings
# valid_column <- function(df, col, no_na = TRUE) {
#   if (!col %in% names(df)) {
#     warning("Column '", col, "' not found in data frame.")
#     return(FALSE)
#   }

#   if (no_na) {
#     invalid <- is.na(df[[col]]) | df[[col]] == ""
#     if (any(invalid)) {
#       return(FALSE)
#       warning(
#         "Column '", col, "' contains NA or empty values at rows: ",
#         paste(which(invalid), collapse = ", ")
#       )
#     } else {
#       return(TRUE)
#     }
#   }

#   TRUE
# }

# Check that a vector column is numeric
is_numeric_col <- function(df, col) {
  if (!col %in% names(df)) {
    warning("Column '", col, "' not found in data frame.")
    return(FALSE)
  }
  if (!is.numeric(df[[col]])) {
    warning("Column '", col, "' is not numeric.")
    return(FALSE)
  }
  TRUE
}

# Check that a list has exactly the expected names
list_has_names <- function(lst, expected_names, warn = TRUE) {
  actual_names <- names(lst)
  missing <- setdiff(expected_names, actual_names)
  extra <- setdiff(actual_names, expected_names)

  if (length(missing) > 0) {
    if (warn) {
      warning("List is missing expected names: ", paste(missing, collapse = ", "))
    }
  }
  if (length(extra) > 0) {
    if (warn) {
      warning("List contains unexpected names: ", paste(extra, collapse = ", "))
    }
  }

  return(length(missing) == 0 && length(extra) == 0)
}

# Check that each element of a list has the expected type
# expected_types should be a named vector/list: names = list names, values = type as string
list_elements_are_types <- function(lst, expected_types, warn = TRUE) {
  all_ok <- TRUE
  for (element_name in names(expected_types)) {
    expected_type <- expected_types[[element_name]]

    if (!inherits(lst[[element_name]], expected_type)) {
      if (warn) warning(
        "Element '", element_name, "' expected type: '", expected_type, "'."
      )
      all_ok <- FALSE
    }
  }

  return(all_ok)
}

# Check that specified columns in a data frame are of expected types
dataframe_cols_are_of_type <- function(df, col_types, na_ok = TRUE) {
  all_ok <- TRUE
  for (col in names(col_types)) {
    expected_type <- col_types[[col]]

    if (!inherits(df[[col]], expected_type)) {
      warning(
        "Column '", col, "' has type '", actual_type,
        "' but expected type '", expected_type, "'."
      )
      all_ok <- FALSE
    }
    if (!na_ok) {
      invalid <- is.na(df[[col]])
      if (any(invalid)) {
        warning(
          "Column '", col, "' contains NA values at rows: ",
          paste(which(invalid), collapse = ", ")
        )
        all_ok <- FALSE
      }
    }
  }
  return(all_ok)
}