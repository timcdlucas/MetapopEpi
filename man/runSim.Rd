% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathogenDiversity.R
\name{runSim}
\alias{runSim}
\title{Run simulation up to certain event number}
\usage{
runSim(pop, end = "end", start = 1)
}
\arguments{
\item{pop}{A MetapopEpi object}

\item{end}{how long to run the simulation until.}

\item{start}{The event number to run the simulation until.}
}
\description{
Each time step, pick a random event and perform that event.
}

