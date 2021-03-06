% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BagHistfp_err.R
\name{BagHistfp_err}
\alias{BagHistfp_err}
\title{Density estimation whit BagHistfp_err}
\usage{
BagHistfp_err(xx, grid, B = 10, dobs)
}
\arguments{
\item{xx}{data vector for Frequency Folygone construction .}

\item{grid}{grid for density evaluation.}

\item{B}{number of Frequency Folygone to aggregate}

\item{dobs}{density values associated whit test sample}
}
\value{
estimation and estimation error.
}
\description{
This function builds "Frequency Folygones"  from the bootstrap samples.
Then averages  the estimates provided  by these estimators and computes their estimation error .
}
