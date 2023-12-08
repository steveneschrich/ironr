#' Extended Distance Matrix Computation
#'
#' @description This function is an extended version of [stats::dist()] which computes
#'   and returns the distance matrix using the specified distance measure to compute
#'   the distances between the rows of a data matrix.
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
dist <- function(
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
  switch(method,
         euclidean = calculate_distance(x, which = which, euclidean_distance, ...),
         pearson = calculate_distance(x, which = which, pearson_distance, ...),
         pearsonm = calculate_distance(x, which = which, pearsonm_distance, ...),
         rmsd = calculate_distance(x, which = which, rmsd_distance, ...),
         presence = calculate_distance(x, which = which, presence_distance, ...),
         multimetric = calculate_distance(x, which = which, multimetric_distance, ...),
         stats::dist(x, method=method, diag=diag,upper=upper,p=p,...)
  )

}

calculate_distance <- function(x, which = "row", ...) {
  if ( which == "row" )
    calculate_distance_by_row(x, ...)
  else
    calculate_distance_by_column(x, ...)
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
calculate_distance_by_row <- function(x, f = euclidean_distance, ...) {
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
    m[i, j] <- m[j, i] <- f(x[i,],x[j, ], ...)
  }
  # Then f for diagonal
  for (i in 1:nrow(x)) {
    m[i, i] <- f(x[i,], x[i,], ...)
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
calculate_distance_by_column <- function(x, f = euclidean_distance, ...) {
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
    m[i, j] <- m[j, i] <- f(x[,i],x[,j], ...)
  }
  # Then f for diagonal
  for (i in 1:ncol(x)) {
    m[i, i] <- f(x[,i], x[,i], ...)
  }

  as.dist(m)
}


#' Use multiple metrics combined as distance metric
#'
#' @details
#' Note that this particular approach can be very misleading, as distance
#' metrics (or measures) can be used and may be on different scales.
#'
#' @param x Vector 1
#' @param y Vector 2
#' @param methods A list of methods (functions) to apply and combine via geometric mean
#'
#' @return The distance between x and y by calculating the geometric mean of
#'  a number of different methods.
#' @export
#'
#' @examples
#'
multimetric_distance <- function(x, y, methods, summarize_by = geometric_mean, ...) {
  stopifnot(all(sapply(methods, is.function)), is.function(summarize_by))

  # Compute all distance metrics
  p <- sapply(methods, function(m) {
    m(x, y, ...)
  })

  summarize_by(p)

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
rmsd_distance <- function(x, y, na.rm = TRUE) {
  stopifnot(is.vector(x),is.vector(y))
  sqrt(mean((x - y)^2, na.rm = na.rm))
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
euclidean_distance <- function(x, y, na.rm = TRUE) {
  stopifnot(is.vector(x),is.vector(y))

  if ( na.rm ) {
    i <- is.na(x) | is.na(y)
    x <- x[!i]
    y <- y[!i]
  }

  sqrt(sum((x-y)^2))
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
pearson_distance <- function(x,y, na.rm = TRUE) {
  stopifnot(is.vector(x),is.vector(y))
  1-cor(x, y, method="pearson", use = "pairwise.complete.obs")
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
#' @param x
#' @param y
pearsonm_distance <- function(x, y, na.rm = TRUE) {
  stopifnot(is.vector(x),is.vector(y))

  r <- pearson_distance(x, y, na.rm = na.rm)
  sqrt(0.5 * r)
}

#' Calculate distance metric based on presence of value
#'
#' @description Calculate the distance between two vectors as the
#' proportion of values non-NA out of the number of non-NA values in
#' the second vector.
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
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
presence_distance <- function(x, y, ...) {
  length(which(!is.na(x) & !is.na(y)))/length(which(!is.na(y)))
}

