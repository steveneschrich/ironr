
#' List CEL files in Directory/Folder
#'
#' This function produces a character vector containing the names of files in the
#' named directory ending in .CEL or .cel.
#'
#' @details
#' @note If the directory does not exist, this function will cause an error indicating
#' this condition (rather than return nothing).
#'
#' @param path The path to list for files
#' @param pattern Any regex pattern to select cel file name (see [stringr::str_detect]).
#' @param ends_with Extension regex (e.g., .CEL)
#' @param ... Any other parameters to pass to [fs::dir_info].
#'
#' @return A character vector of filenames
#' @export
#'
#' @examples
list.celfiles <- function(path = ".", pattern = NULL, ends_with = "\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$",...) {

  files <- fs::dir_info(path = path, type = c("file","symlink"), regexp = ends_with, ...)

  if ( !is.null(pattern) ) {
    files <- dplyr::mutate(files, file = basename(path))
    files <- dplyr::filter(files, stringr::str_detect(file, pattern))
  }

  files$path
}

#' Read CEL files into an AffyBatch
#'
#' @param dir Find CEL files in dir to read
#' @param pattern Any CEL file pattern (regex): see [stringr::str_detect]
#' @param ... Any other parameters to [affy::read.affybatch]
#'
#' @return An [affy::AffyBatch] object.
#' @export
#'
#' @examples
read.affybatch <- function(dir, pattern = NULL,...) {
  suppressWarnings({
    affy::read.affybatch(filenames = list.celfiles(dir, pattern=pattern), ...)
  })
}

#' Read CEL files into an Oligo ExpressionFeatureSet
#'
#' @param dir Find CEL files in dir to read
#' @param pattern Any CEL file pattern (regex): see [stringr::str_detect]
#' @param ... Any other parameters to [oligo::read.celfiles]
#'
#' @return An [oligoClasses::ExpressionFeatureSet]
#' @export
#'
#' @examples
read.oligo <- function(dir, pattern = NULL, ...) {
  oligo::read.celfiles(filenames = list.celfiles(dir, pattern=pattern), ..., verbose = FALSE)
}


#' Extract perfect match (PM) probes
#'
#' This function extracts PM probes from the object based on the object type.
#'
#' @param x An AffyBatch or ExpressionFeatureSet
#'
#' @return A matrix of PM values
#' @export
#'
#' @examples
perfect_match_probes <- function(x) {
  if (methods::is(x, "AffyBatch"))
    return(affy::pm(x))
  if (methods::is(x, "FeatureSet"))
    return(oligo::pm(x))
}
