# NB: This would be the driver function
# It would need to load in the cel files (optionally affybatch, or just for
# the moment pm array). Then log transform, remove low values to NA (optionally),
# then calculate using rmsd (or other options).
#
# This almost works. The values are off by a bit, but generally it is close.

iron_findmedoid <- function() {

  cf <- list.files(path=here::here("tests"),pattern = "^T.*\\.CEL", full.names = TRUE)

  # Not a fan of suppressing warnings, but the phenodata it creates takes a bit
  # of effort, so this is easier to maintain compatability with Bioconductor.
  ab <- suppressWarnings(
    affy::read.affybatch(
      filenames = cf,
      #phenoData = Biobase::AnnotatedDataFrame(data.frame(files=cf)),
      rm.mask = FALSE, # libaffy setting
      rm.outliers = FALSE, # libaffy setting
      verbose = TRUE,
      sd = FALSE # sd is not needed
    )
  )

  # Note, right now we use libaffy::get_pm() for the rm.control, but ultimately
  # the iron version of find median will remove it.
  pm_array <- log2(libaffy::get_pm(ab, rm.control = TRUE))
  pw_distances <- pairwise_col_distance(pm_array, distance_method = "rmsd",overlap_normalization=TRUE)

  d <- data.frame(
    median_sample = ifelse(pw_distances == min(pw_distances), "*",""),
    sample_name = names(pw_distances),
    mean_distance = pw_distances,
    mean_intensity = colMeans(pm_array),
    row.names = NULL
  )


  print(d)
  cat("Average Distance:" , mean(pw_distances),"\n")

  invisible(d)
}



# These have to be totally,super hacky versions to accomodate the needs
# Also, we'll need a Rscript/getopt wrapper as well. These should
# probably just be moved over to "findmedian.R" which is what it is.
findmedian.matrix <- function(x, ignore_weak = TRUE, distance = "rmsd", verbose = FALSE, ...) {

  # Calculate distances between columns using standard approach. Note that we
  # do not use a adjustment to the distance for missingness, since this added later
  # (asymmetrically)
  d <- distance(x, which = "col", method = distance, na.action = NULL, ...)

  # Convert dist to a matrix. The distances must be symmetric unless the distance
  # function returns a full matrix (not a dist object). Currently it does not. We
  # add asymmetry below.
  d <- as.matrix(d)

  # NA the diagonals. This is not needed for minimizing the sum, but mean values
  # should not include the distance between self (even if it's 0).
  diag(d) <- NA

  # Use custom distance normalization for missingness (asymmetric, see the function).
  offdiag <- utils::combn(1:ncol(x), 2, simplify = FALSE)
  for (pair in offdiag ) {
    i <- pair[1]
    j <- pair[2]
    d[i, j] <- d[i, j] * na.dist_norm_iron(x[,i], x[,j])
    d[j, i] <- d[j, i] * na.dist_norm_iron(x[,j], x[,i])
  }


  # Calculate the sums (medoid definition) and find the minimum sum of
  # distances. Note that we use the rowSums in either case, under the assumption
  # that if there is asymmetry the row-wise distances are appropriate.
  dist.summary <- rowMeans(d, na.rm = TRUE)
  mi <- which.min(dist.summary)

  # Finally, the result is either the row or column in question.
  e <- exemplar(x, index = mi, which = "col", dist.summary = dist.summary,
           dist.method = distance, summary.method = "mean")

  if ( verbose )
    print(e)

  e

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

findmedian.affybatch <- function() {libaffy::getpm()}






