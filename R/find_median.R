#' Find median sample from dataset
#'
#' Implements a similar function to libaffy IRON findmedian function
#'
#' @details
#' This function approximates the findmedian program from libaffy. Briefly,
#' a medoid is identified using the rmsd distance metric and finding the
#' sample with the smallest all-sample mean distance. There are numerous
#' differences with the C code which is highly specific, but this function
#' roughly approximates the functionality.
#'
#' Note that a more complete view of this process is given in the vignette
#' or in the functions [medoid()] or [find_exemplar()] in this package.
#'
#'
#' @param x
#' @param which
#' @param na.rm
#' @param distance
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
find_median <- function(x, which=c("col","row"), na.rm = TRUE, distance = "rmsd", ...) {
  find_exemplar(x, which=which, distance=distance, summary.method=mean)
}

#' Calculate conditional probability of overlap given query
#'
#' @description Calculate the conditional probability between two vectors as the
#' proportion non-NA values in the second vector by the total number of non-NA values
#' in common.
#'
#' @details
#'
#' Conditional probability of #overlap given #query.
#'
#' This is asymmetric, but that's what we want in this case.
#' We are searching for candidate reference samples.
#' Say that we have a reference with 100% present, and the query
#' is only 25% present.  We don't want to use some similarity metric
#' that combines those two #present together, since we want to know
#' if the query is maximally covered or not.  #overlap / #query will
#' be 1.0, since it is fully covered.  The other way around,
#' #overlap / #target, would be 0.25.  This is the scoring behavior
#' we want, since the first case indicates a good reference sample
#' and the 2nd case indicates a poorer reference sample.
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
na.dist_norm_iron <- function(x,y) {
  present_in_both <- !(is.na(x) | is.na(y))
  present_in_y <- !is.na(y)


  sum(present_in_y) / sum(present_in_both)
}







