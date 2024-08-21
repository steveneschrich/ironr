#' Title
#'
#' @param dir
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
iron.cel.findmedian <- function(dir, pattern=NULL) {

  x <- read.oligo(dir, pattern=pattern)
  xm <- log2(perfect_match_probes(x))
  xm[xm<=0] <- NA
  ex <- find_median(xm)

  ex
}

#' Title
#'
#' @param dir
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
iron.cel.medoid<- function(dir, pattern=NULL) {
  x <- read.oligo(dir, pattern=pattern)
  xm <- log2(perfect_match_probes(x))
  xm[xm<=0]<-NA
  ex <- medoid(xm)

  ex
}
