#' Extended Distance Matrix Computation
#'
#' @description This function is an extended version of [stats::dist()] which computes
#'   and returns the distance matrix using the specified distance measure to compute
#'   the distances between the rows or columns of a data matrix.
#' @param x
#' @param method
#' @param diag
#' @param upper
#' @param p
#'
#' @return
#' @export
#'
#' @examples
distance <- function(
    x,
    method = c("euclidean","pearson","pearsonm","rmsd","multimetric","maximum","manhattan","canberra",
               "binary","minkowski", "presence"
    ),
    diag = FALSE,
    upper = FALSE,
    p = 2,
    which = c("row","column"),
    ...
) {

  method <- match.arg(method)
  which <- match.arg(which)

  if ( method %in% c("maximum","manhattan","canberra",
                     "binary","minkowski") && which == "column") {
    stop("Default distance methods not implemented via column yet.")
  }

  if ( which == "row" ) {
    d <- switch(method,
           euclidean = rowDistance(x, euclidean_distance),
           pearson = rowDistance(x,  pearson_distance),
           pearsonm = rowDistance(x, pearsonm_distance),
           rmsd = rowDistance(x,  rmsd_distance),
           stats::dist(x, method=method, diag=diag,upper=upper,p=p)
    )
  } else {
    d <- switch(method,
           euclidean = colDistance(x, euclidean_distance),
           pearson = colDistance(x,  pearson_distance),
           pearsonm = colDistance(x,  pearsonm_distance),
           rmsd = colDistance(x, rmsd_distance),
           stats::dist(x, method=method, diag=diag,upper=upper,p=p)
    )
  }

  d
}

#' Calculate distances between rows in a matrix
#'
#' @description Calculate distances between rows of a matrix using the given
#' distance metric (as a function) and return a dist object.
#'
#' @param x Matrix with multiple rows
#' @param f A function to use for calculating distances
#' @param ... Additional parameters needed for f
#'
#' @return A dist object representing distances between rows.
#' @export
#'
#' @examples
rowDistance <- function(x, f = euclidean_distance) {
  x <- as.matrix(x)
  # Initialize array of results
  m <- matrix(
    data = 0,
    nrow = nrow(x), ncol = nrow(x),
    dimnames = list(rownames(x), rownames(x))
  )
  # Calculate pairs for computing f
  offdiag <- utils::combn(1:nrow(x), 2, simplify = FALSE)
  # Call f for all pairs
  for (pair in offdiag) {
    i <- pair[1]
    j <- pair[2]
    m[i, j] <- m[j, i] <- f(x[i,],x[j, ])
  }
  # Then f for diagonal
  for (i in 1:nrow(x)) {
    m[i, i] <- f(x[i,], x[i,])
  }

  as.dist(m)
}



#' Calculate distances between columns in a matrix
#'
#' @description Calculate distances between columns of a matrix using the given
#' distance metric (as a function) and return a dist object.
#'
#' @param x Matrix with multiple columns
#' @param f A function to use for calculating distances
#' @param ... Additional parameters needed for f
#'
#' @return A dist object representing distances between columns.
#' @export
#'
#' @examples
colDistance <- function(x, f = euclidean_distance) {
  x <- as.matrix(x)
  # Initialize array of results
  m <- matrix(
    data = 0,
    nrow = ncol(x), ncol = ncol(x),
    dimnames = list(colnames(x), colnames(x))
  )
  # Calculate pairs for computing f
  offdiag <- utils::combn(1:ncol(x), 2, simplify = FALSE)
  # Call f for all pairs
  for (pair in offdiag) {
    i <- pair[1]
    j <- pair[2]
    m[i, j] <- m[j, i] <- f(x[,i],x[,j])
  }
  # Then f for diagonal
  for (i in 1:ncol(x)) {
    m[i, i] <- f(x[,i], x[,i])
  }

  as.dist(m)
}




#' Calculate rmsd between two vectors
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
rmsd_distance <- function(x, y) {
  stopifnot(is.vector(x),is.vector(y))


  d <- mean((x - y)^2, na.rm = na.rm)

  d <- d * na.dist_norm(x,y)

  sqrt(d)
}

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
euclidean_distance <- function(x, y) {
  stopifnot(is.vector(x),is.vector(y))

  i <- is.na(x) | is.na(y)

  d <- sum(( x[!i] - y[!i] )^2)

  d <- d * na.dist_norm(x,y)

  sqrt(d)
}

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
pearson_distance <- function(x,y) {
  stopifnot(is.vector(x),is.vector(y))
  d <- 1-cor(x, y, method="pearson", use = "pairwise.complete.obs")


    d <- d * na.dist_norm(x, y)

  d
}

#' Calculate distance metric based on Pearson correlation
#'
#' @description Calculate distance between two vectors using 1-Pearson correlation
#' with a correction to make a metric.
#'
#' @details
#'
#'
#' transform r into a metric distance [sqrt(0.5 * (1-r))] */
#' see Stijn van Dongen and Anton J. Enright 2012
#' "Metric distance derived from cosine similarity and Pearson
#'      *  and Spearman correlations"
#'
#' @param x Vector
#' @param y Vector
pearsonm_distance <- function(x, y) {
  stopifnot(is.vector(x),is.vector(y))

  r <- pearson_distance(x, y)
  sqrt(0.5 * r)
}




#' Calculate distance normalization factor based on missingness
#'
#' @description Calculates the distance normalization factor (NORMAL) for
#'  missing data (see details). This is a multiplicative factor which
#'  when multiplied by distance, increases distances by the number of
#'  missing values.
#'
#' @details When calculating distances between two vectors, the number of
#' missing values in these vectors can impact on the distance. With too few
#' non-missing values, the distance may be inflated.
#'
#' This function is a helper function that be used in the [distance()] functions
#' in this package to include a distance penalty for missingness. Specifically,
#' each distance function provides a `na.action` parameter which will multiply
#' a distance penalty (as implemented here) with the distance calculation.
#'
#' This function implements a common approach that is given in:
#' \preformatted{
#' John K. Dixon
#' "Pattern Recognition with Partly Missing Data"
#' IEEE Transactions on Systems, Man, and Cybernetics
#' Volume: 9, Issue: 10, pp. 617 - 621, Oct. 1979.
#' http://ieeexplore.ieee.org/abstract/document/4310090/
#' }
#'
#' In this approach (defined for squared Euclidean distance), the normalization
#' factor is
#' \eqn{\dfrac{N}{N-NB}}
#' where `N` is the total number of observations and `NB` is the number of blanks.
#'
#' Note that when computing the Euclidean distance (as was done in the paper),
#' the normalization factor is multiplied by the sum of squared differences. Then
#' the square root could be applied to the entire term. The assumption for this
#' function is that any modification is d
na.dist_norm <- function(x, y) {
  i <- is.na(x) | is.na(y)
  length(x) / sum(!i)
}


