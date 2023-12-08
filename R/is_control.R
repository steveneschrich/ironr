#' Determine if probe name is control
#'
#' @description Determine which of a vector of probe names is a control probe
#'  based on heuristics from the probe name.
#'
#' @details
#' Affymetrix probes/probesets have a variety of meanings embedded in the
#' names. In this case, there are certain probes on the array that translate
#' to control-related activities. These include:
#'   - Starts with AFFX
#'   - Includes the text 'spikein' in a variety of formats
#'   - Includes the word 'control'
#'
#' This function detects these various heuristics and indicates which of the
#' names are control probes.
#'
#' @param p A character vector of probe names
#'
#' @return Logical value indicating if the name(s) are controls
#' @export
#'
#' @examples
#' affy_is_control(c("AFFX-12345","12345_at","This_is_a_spike-in"))
#'
affy_is_control <- function(p) {
  # Imported logic from libaffy:
  # https://gitlab.moffitt.usf.edu:8000/WelshEA/libaffy/-/blob/master/libaffy/utils/is_control_probe.c

  # Starts with AFFX
  stringi::stri_startswith_fixed(p, "AFFX", case_insensitive = TRUE) |

    # Includes spike[-_ ]in (or no space)
    stringi::stri_detect_regex(p, "spike[-_ ]?in") |

    # Includes term control
    stringi::stri_detect_fixed(p, "control", case_insensitive = TRUE)
}
