#NB: There should be a log2 check.

#
# NB: How are use_exclusions and use_spikeins used?


find_median.matrix<-function(X, positive_only=TRUE,
                             mean_center=TRUE,
                             use_log=FALSE,
                             min_threshold=1,
                             ignore_weak=FALSE) {

  # Start manipulating the input relative to flags
  Xp<-X

  # ignore_weak means set values < 0 to NA
  if (ignore_weak) {
    message("Ignoring points with weak values.")
    Xp[which(Xp<0)]<-NA
  }

  # logs are ok, after truncating at min_threshold
  if (use_log) {
    message("Log2 transforming data.")
    Xp[which(Xp<min_threshold)]<-min_threshold
    Xp<-log2(Xp)
  }

  # Mean-center columns
  if (mean_center) {
    message("Mean-centering vectors.")
    Xp<-scale(Xp, scale=FALSE, center=TRUE)
  }


  # Q for Eric: do we remove negatives after scaling??
  if ( positive_only) {
    X[which(X<0)]<-NA
  }

  return(pairgen_find_median_chip_distance(Xp))
}
