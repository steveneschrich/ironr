# normalize.R
# normalization and scaling routines

normalize_mean_center <- function(x, na.rm = TRUE, trim = 0) {
  x <- as.matrix(x)

  if ( trim != 0 ) {
    xp <- apply(x, 2, \(.x) { mean(.x, na.rm = na.rm, trim = trim)})
  } else {
    xp <- colMeans(x)
  }

  scale(xp, center = xp, scale = FALSE)

}
