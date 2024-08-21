

exemplar <- function(x = NULL, index = 1, which = c("row","col"),
                     dist.summary = NULL, dist.method = "", summary.method="") {
  if ( is.null(x) ) {
    stop("x must be a data frame or matrix of data.")
  }
  new_exemplar(x, index, which, dist.summary, dist.method, summary.method)
}

new_exemplar <- function(x = NULL, index, which = c("row","col"),
                     dist.summary = NULL, dist.method = "", summary.method) {
  which <- match.arg(which)
  # Error checking and extract exemplar
  stopifnot(is.matrix(x) || is.data.frame(x))
  if ( which == "row" ) {
    stopifnot(index > 0, index <= nrow(x))
    e <- x[index,]
  } else {
    stopifnot(index > 0, index <= ncol(x))
    e <- x[,index]
  }
  # Create the structure
  e <- structure(
    e,
    x = x,
    index = index,
    which = which,
    dist.summary = dist.summary,
    dist.method =  dist.method,
    summary.method = summary.method,
    class = "exemplar"
  )


  e
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.exemplar <- function(x, ...) {
  print(format.exemplar(x, ...))
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.exemplar <- function(x, ...) {
  print(names(attr(x, "index")))
}
#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
format.exemplar <- function(x, ...) {

  which <- attr(x, "which")
  data <- attr(x, "x")

  if ( which == "col") {
    o <- data.frame(
      medoid = rep("",ncol(data)),
      names = colnames(data),
      column_mean = colMeans(data, na.rm = TRUE),
      distance = attr(x, "dist.summary")
    )
  } else {
    o <- data.frame(
      medoid = rep("",nrow(data)),
      names = rownames(data),
      row_mean = rowMeans(data, na.rm = TRUE),
      distance = attr(x, "dist.summary")
    )
  }
  o[attr(x, "index"), "medoid"] <- "*"

  rownames(o) <- NULL

  o
}



#' Find exemplar of a matrix
#'
#' @param x Matrix to operate on
#' @param which One of c("cols","rows") for the exemplar
#' @param distance A distance measure to use (see [distance()] for details).
#' @param summary.method A summary method applied to distances (default is sum)
#'  which is then minimized to find the exemplar.
#' @param ... Any further parameters to distance function.
#'
#' @description Identify the exemplar (e.g., a medoid) of the matrix.
#'
#' @details This function is a generic approach to identifying a medoid or
#'  exemplar in a matrix. A medoid is the column (row) with the minimum
#'  summed distance to all other columns (rows). However, it is possible to
#'  consider other operations than sum of distances (e.g., mean distance). This
#'  function provides the non-standard implementation of medoid (exemplars).
#'
#'  While this function is only slightly more complex, see [medoid()] for
#'  the definition. The implementation is slightly more efficient than this one,
#'  due to the generic nature of the summary method.
#'
#'  Note that one of the motivating factors for this extension to medoids is
#'  that the IRON algorithm uses the minimum **mean** distance vs. the minimum
#'  **summed** distance.
#'
#' @export
find_exemplar <- function(x, which = c("col","row"),
                                distance = "euclidean",
                                summary.method = sum,
                                ...) {
  which <- match.arg(which)

  # Calculate the distance between columns/rows as a 'dist' object
  d <- distance(x, which = which,  na.rm = na.rm, method = distance, ...)

  # Convert dist to a matrix. The distances must be symmetric unless the distance
  # function returns a full matrix (not a dist object). Currently it does not.
  d <- as.matrix(d)

  # NA the diagonals. This is not needed for minimizing the sum, but other summary
  # functions should not include the distance between self (even if it's 0).
  diag(d) <- NA


  # Calculate the sums (medoid definition) and find the minimum sum of
  # distances. Note that we use the rowSums in either case, under the assumption
  # that if there is asymmetry the row-wise distances are appropriate.
  dist.summary <- apply(d, 1, \(.x) {
    .x <- na.omit(.x)
    summary.method(.x)
  })
  mi <- which.min(dist.summary)



  # Finally, the result is either the row or column in question.
  exemplar(x, index = mi, which = which, dist.summary = dist.summary,
           dist.method = distance,summary.method = summary.method)
}

#' Find the medoid of a matrix
#'
#' @description Find the column (row) of the matrix with the minimum summed distance
#'  to all other columns (rows).
#'
#' @details The medoid is defined as the example with the minimum summed distances to
#' all other examples in the dataset. This function calculates all pair-wise distances
#' between columns (or rows) and then identifies the columns (or row) with the minimum
#' distance.
#'
#' Some notes:
#'
#' - If multiple columns (rows) all have a minimum summed distance, then the first column (row) with
#' the identical minimum distance is selected (see [base::which.min()]).
#'
#' - The distance measure must be specified as a string. See [distance()] in this package
#' for details on available options.
#'
#' - An assumption is made in this function that the distance measure is symmetric. That is,
#' the distance between `A` and `B` is the distance between `B` and `A`.
#'
#' @param x A matrix to compare columns or rows
#' @param which Identify the medoid "row" or "col"
#' @param na.rm Should NA be removed prior to distance calculation (default TRUE)
#' @param distance A string representing a distance measure to use (see [distance()])
#' @param ... Any parameters to pass to the distance function
#'
#' @return An object of type [exemplar()] representing the exemplar.
#' @export
#'
#' @examples
medoid <- function(x, which = c("col","row"), na.rm = TRUE, distance = "euclidean", ...) {
  which <- match.arg(which)

  # Calculate the distance between columns/rows as a 'dist' object
  d <- distance(x, which = which,  na.rm = na.rm, method = distance, ...)

  # Convert dist to a matrix. The distances must be symmetric unless the distance
  # function returns a full matrix (not a dist object). Currently it does not.
  d <- as.matrix(d)

  # Calculate the sums (medoid definition) and find the minimum sum of
  # distances. Note that we use the rowSums in either case, under the assumption
  # that if there is asymmetry the row-wise distances are appropriate.
  dist.summary <- rowSums(d, na.rm = na.rm)
  mi <- which.min(dist.summary)

  # Finally, the result is either the row or column in question.
  exemplar(x, index = mi, which = which, dist.summary = dist.summary,
           dist.method = distance,summary.method = "sum")

}

