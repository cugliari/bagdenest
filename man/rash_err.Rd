% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rash_err.R
\name{rash_err}
\alias{rash_err}
\title{density estimation and estimation error whith rash.}
\usage{
rash_err(xx, grid, nbr = 50, B = 10, dobs)
}
\arguments{
\item{xx}{data vector}

\item{grid}{grid for density evaluation}

\item{nbr}{number of breaks for histogram}

\item{B}{number of histograms to aggregate}

\item{dobs}{density values associated whit test sample}
}
\value{
estimation and estimattion error
}
\description{
this function builds several histograms in disturbing by translation an initial partition
             of the histogram with a Gaussian noise and then aggregates them.
}
