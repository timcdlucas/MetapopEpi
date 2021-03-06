% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathogenDiversity.R
\name{waitingTime}
\alias{waitingTime}
\title{Takes a number [0,1] (which is presumably random) and returns waiting time.}
\usage{
waitingTime(pop, t)
}
\arguments{
\item{pop}{A population object}

\item{t}{Time step. Note the population should be AT time t, going to t+1.}
}
\description{
A single draw from a exponential distribution given randU, a number between 0 and 1.
  As population now has
}

