tukey_biweight <- function(x, c = 5, epsilon = 1e-04) {
  if ( !is.vector(x) ) {
    cli::cli_abort("x is not a vector, cannot calculate tukey biweight.")
  }
  # X-M (need absolute and not)
  x_minus_m <- x - median(x, na.rm = TRUE)
  # Calculate S, median |X-M|
  S <- median(abs(x_minus_m))
  # Calculate distances u
  u <- x_minus_m/(c * S + epsilon)

  # w(u) = (1-u^2)^2 unless |u| > 1 => 0
  # We just replace |u| > 1 with 1 to get 0's from the following
  u[abs(u)>1] <- 1
  w <- (1 - (u^2))^2
  # Tbi = sum(w*x)/sum(w)
  T_bi <- sum(w*x)/sum(w)

  return(T_bi)
}


#' Scaled-down calculation of probeset using MAS5
#'
#' @description Internal IRON function for calculating MAS5 probeset signal
#' without additional processing (following libaffy).
#'
#' @details
#'
#' @note This function operates on the probeset level, meaning the input
#' should either be a vector of probe values for a probeset or a matrix
#' of probe values (row) for a set of samples (cols).
#'
affy_mas5_signal <- function(cels, delta = 2.0e-20, keep_log2 = FALSE) {
  probe_mapping <- affy::indexProbes(cels, which="pm")

  # Extract all intensities
  celdata <- affy::intensity(cels)

  # Floor data and log2 transform
  x[x<delta]<-delta
  x <- log2(x)


  # Iterate over samples (columns), then calculate summary for each probeset.
  signal <- apply(celdata, 2, function(celfile) {
    sapply(
      vctrs::vec_chop(celfile, indices = probe_mapping),
      tukey_biweight
    )
  })
  # The probeset names aren't automatically included
  rownames(signal) <- names(probe_mapping)

  # Make sure the result is a matrix
  signal <- as.matrix(signal)

  if ( !keep_log2 )
    signal <- 2^signal

  signal
}

