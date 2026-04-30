#' A sample `data.frame` containing Differential Expression Analysis, generated
#' with `limma`
#'
#' @details This `data.frame` object contains the results of a Differential
#' Expression Analysis performed (with `limma`) on data from the `macrophage`
#' package, contrasting the counts from naive macrophage to those associated
#' with IFNg.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the KEGGemUP package, the file is called `create_deresults.R`.
#'
#' @return A sample `data.frame` object, extracted from the `topTable` function
#' as a result of running the `limma` DE workflow
#'
#' @format A `data.frame` object
#'
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name res_de_macro_IFNg_vs_naive
#' @docType data
NULL
