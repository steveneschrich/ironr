% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance.R
\name{rowDistance}
\alias{rowDistance}
\title{Calculate distances between rows in a matrix}
\usage{
rowDistance(x, f = euclidean_distance)
}
\arguments{
\item{x}{Matrix with multiple rows}

\item{f}{A function to use for calculating distances}

\item{...}{Additional parameters needed for f}
}
\value{
A dist object representing distances between rows.
}
\description{
Calculate distances between rows of a matrix using the given
distance metric (as a function) and return a dist object.
}
