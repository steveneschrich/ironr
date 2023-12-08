
# NB: This would be the driver function
# It would need to load in the cel files (optionally affybatch, or just for
# the moment pm array). Then log transform, remove low values to NA (optionally),
# then calculate using rmsd (or other options).
#
# This almost works. The values are off by a bit, but generally it is close.

iron_findmedian <- function() {
  pm_array <- affyio::read.celfile.probeintensity.matrices(
    filenames=list.files(path=here::here("tests"),pattern = "^T.*\\.CEL", full.names = TRUE),
    cdfInfo = as.list(hgu133plus2cdf::hgu133plus2cdf),
    verbose = TRUE
  )$pm

  pm_array <- log2(pm_array)
  d <- median_col_distance(pm_array, distance_method = "rmsd")
  cat(d)
  cat("Average RMSD:", mean(d$mean_distance))
}


# this is actually (at it's heart) a bit simpler. In terms of matrices, etc
# it is to calculate column distances (whoops, I did row distances elsewhere).
# Column distances, and then find the median one somehow. right?
# median_col_distance(verbose = TRUE)
#.  distance matrix
#.  median (as pair?)
# Conditional probability of #overlap given #query.
#
# This is asymmetric, but that's what we want in this case.
# We are searching for candidate reference samples.
# Say that we have a reference with 100% present, and the query
# is only 25% present.  We don't want to use some similarity metric
# that combines those two #present together, since we want to know
# if the query is maximally covered or not.  #overlap / #query will
# be 1.0, since it is fully covered.  The other way around,
# #overlap / #target, would be 0.25.  This is the scoring behavior
# we want, since the first case indicates a good reference sample
# and the 2nd case indicates a poorer reference sample.
#

median_col_distance <- function(x, verbose = TRUE, na.rm = TRUE,
                                distance_method = "euclidean",
                                overlap_normalization = FALSE,
                                ...) {

  col_distances <- dist(x, which = "col", na.rm =TRUE, method = distance_method, ...)

  # This is a particular variation which normalizes distances by the total
  # amount of overlap.
  if ( overlap_normalization ) {
    col_distances <- col_distances / dist(which="col", method = "presence")
  }

  # Return a summary of distances
  d <- data.frame(
    name = colnames(x),
    mean_distance = apply(as.matrix(col_distances),2,mean),
    column_mean = apply(x, 2, \(.x) { mean(.x, na.rm = TRUE)}),
    row.names = NULL
  )

  d

}


pairgen_find_median_chip_distance<-function(X) {
  nc <- ncol(X)
  x_names<-colnames(X)

  # Initialize empty matrices for each result.
  rmsd <- pearson <- geom_mean <- matrix(0, nrow=nc, ncol=nc,
                                         dimnames=list(x_names, x_names))
  # This is how cor does it, calculate the lower triangle.
  for (i in seq_len(nc) ) {
    for (j in seq_len(i)) {
      rmsd[i, j] <- calculate_rmsd(X[,i], X[,j])
      pearson[i, j] <- calculate_pearson_distance(X[,i], X[,j])
      geom_mean[i, j] <- sqrt(rmsd[i, j] * pearson[i, j])
      #normalization_factor[i, j] <- length(complete.cases(X[,i], X[,j]))
    }
  }

  # NB: I replace the diagonal with NA because I don't want to consider
  # self-similarity here.
  diag(rmsd)<-NA
  diag(pearson)<-NA
  diag(geom_mean)<-NA

  # Normalize distances by pmissing
  #/* Conditional probability of #overlap given #query.
  #*
  #  * This is assymetric, but that's what we want in this case.
  #       * We are searching for candidate reference samples.
  #       * Say that we have a reference with 100% present, and the query
  #       * is only 25% present.  We don't want to use some similarity metric
  #* that combines those two #present together, since we want to know
  #* if the query is maximally covered or not.  #overlap / #query will
  # be 1.0, since it is fully covered.  The other way around,
  #* #overlap / #target, would be 0.25.  This is the scoring behavior
  #  * we want, since the first case indicates a good reference sample
  #* and the 2nd case indicates a poorer reference sample.
  #*/
  # stuff<-normalization_factor (after filling it out) / whichever
  #overlap/#query
  # complete.cases between a,b
  #


  # The lower triangle is calculated, reflect that to the upper with
  # a transpose. This is what cor() does.
  # This means that adding the transpose doesn't overcount, since
  # the NA's are propogated.
  rmsd <- rmsd + t(rmsd)
  pearson <- pearson + t(pearson)
  geom_mean <- geom_mean + t(geom_mean)

  # Calculate mean distance to all other chips
  rmsd_mean<-apply(rmsd,1, mean, na.rm=T)
  pearson_mean<-apply(pearson, 1, mean, na.rm=T)
  geom_mean<-apply(geom_mean, 1, mean, na.rm=T)

  # For now, print out the same output as Eric.
  print(data.frame(Score=rep("Score",length(rmsd_mean)),
                   SampleIndex=0:(length(rmsd_mean)-1),
                   MeanDistance=rmsd_mean,
                   SampleName=names(rmsd_mean),
        MeanLog2Abundance=colMeans(X,na.rm=TRUE)))
  cat(sprintf("Average RMSD: %f\n",mean(rmsd_mean)))
  # Return the mean distances.

}


