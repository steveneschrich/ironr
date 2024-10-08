% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance.R
\name{na.dist_norm}
\alias{na.dist_norm}
\title{Calculate distance normalization factor based on missingness}
\usage{
na.dist_norm(x, y)
}
\description{
Calculates the distance normalization factor (NORMAL) for
 missing data (see details). This is a multiplicative factor which
 when multiplied by distance, increases distances by the number of
 missing values.
}
\details{
When calculating distances between two vectors, the number of
missing values in these vectors can impact on the distance. With too few
non-missing values, the distance may be inflated.

This function is a helper function that be used in the [distance()] functions
in this package to include a distance penalty for missingness. Specifically,
each distance function provides a `na.action` parameter which will multiply
a distance penalty (as implemented here) with the distance calculation.

This function implements a common approach that is given in:
\preformatted{
John K. Dixon
"Pattern Recognition with Partly Missing Data"
IEEE Transactions on Systems, Man, and Cybernetics
Volume: 9, Issue: 10, pp. 617 - 621, Oct. 1979.
http://ieeexplore.ieee.org/abstract/document/4310090/
}

In this approach (defined for squared Euclidean distance), the normalization
factor is
\eqn{\dfrac{N}{N-NB}}
where `N` is the total number of observations and `NB` is the number of blanks.

Note that when computing the Euclidean distance (as was done in the paper),
the normalization factor is multiplied by the sum of squared differences. Then
the square root could be applied to the entire term. The assumption for this
function is that any modification is d
}
