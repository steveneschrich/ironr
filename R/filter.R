#' Filter observations based on signal
#'
#' @description Filter observations in two vectors based on IRON definitions.
#'
#' @details
#' This filter implements logic associated with the IRON approach to normalization.
#' The logic is not completely transparent, but is included in text below:
#'
#' - Exclude missing values in x or y
#' - Exclude values <= min_value
#' - Exclude values >= max_value
#' - Exclude values that are the minimum (>= min_value) of the vector
#'
#' If these conditions result in no values being available (all filtered out),
#' then only the top 2 conditions are used for filtering.
#'
#'
#'
#' @param x Vector
#' @param y Vector
#' @param min_value Filter out values <= min_value (default 1e-05)
#' @param max_value Filter out values >= max_value (default Inf)
#'
#' @return A vector of logical values indicating whether or not *keep* values
#'   based on parallel comparison of x and y.
#' @export
#'
#' @examples
iron_filter2 <- function(x, y, min_value = 1e-05, max_value = Inf, verbose = TRUE) {

  # Remove values below the min_value
  x[x < min_value] <- NA
  y[y < min_value] <- NA

  # These are the base conditions to exclude:

  # Missing values
  exclude_missing <- is.na(x) | is.na(y)
  # Values at minimum value (below minimum already handled above)
  exclude_min <- !exclude_missing & (( x == min_value ) | (y == min_value))
  # Values above maximum value
  exclude_max <- !exclude_missing & (
    ( x >= max_value) | ( y >= max_value)
  )
  # This one is a little strange, exclude the minimum value that is at or above the min_value.
  # Since we have NA'd values less than min_value, whatever left (possibly min_value) is the
  # minimum.
  if ( all(exclude_missing ) ) {
    exclude_min_observed_value <- exclude_missing
  } else {
    exclude_min_observed_value <- !exclude_missing & (
      x <= min(x, na.rm = TRUE) |
        y <= min(y, na.rm = TRUE)
    )
  }
  # The purpose of piecewise filtering is to only filter on max/min_observed_value if we
  # still have values to filter on.
  exclusions <- ( exclude_missing | exclude_min | exclude_max | exclude_min_observed_value )
  if ( all(exclusions) ) {
    cli::cli_alert_warning("Using less strict filtering criteria to avoid empty set.")
    exclusions <- ( exclude_missing | exclude_min )
  }

  !exclusions

}


