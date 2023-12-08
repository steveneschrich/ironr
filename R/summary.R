#' Calculate geometric mean of vector
#'
#' @param x Vector
#' @param na.rm (TRUE) Remove NA before calculating
#' @param positive_only (TRUE) Are only positive values used in geometric mean?
#'
#' @return The geometric mean of the input vector
#' @export
#'
#' @examples
geometric_mean <- function(x, na.rm = TRUE, positive_only=TRUE) {
  if ( positive_only ) {
    x[x <= 0] <- NA
  }
  exp(mean(log(x), na.rm = na.rm))
}
