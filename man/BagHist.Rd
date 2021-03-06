% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BagHist.R
\name{BagHist}
\alias{BagHist}
\title{Density estimation whit BagHist}
\usage{
BagHist(xx, grid, B = 10)
}
\arguments{
\item{xx}{data vector for histograms construction .}

\item{grid}{grid for density evaluation}

\item{B}{number of histograms to aggregate}
}
\value{
estimation
}
\description{
This function builds histograms from the bootstrap samples.
Then averages  the estimates provided  by these histograms and computes thes estimation errro.
}
