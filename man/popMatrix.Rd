% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathogenDiversity.R
\name{popMatrix}
\alias{popMatrix}
\title{Calculate a dispersal matrix}
\usage{
popMatrix(locations, maxDistance, kern)
}
\arguments{
\item{locations}{x, y numeric locations of each colony.}

\item{maxDistance}{For kernels with an upper threshold}

\item{kern}{Character, giving the name of the distance kernel to be used.}
}
\description{
Calulate a weighted (or unweighted) adjacency matrix. Then weight this value to convert to a dispersal matrix.
  An element i, j is the proportion of individuals dispersing from i that go to j.
}

