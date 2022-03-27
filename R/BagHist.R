#' @title Density Estimation with Bagged Histograms (BagHist)
#'
#' @description Estimate the density of a univariate data vector using BagHist. First,
#' obtain bootstrap samples of the input vector data. For each sample, a histogram is
#' computed using an optimized number of break points. The BagHist estimation is the mean of the individual histograms.
#'
#' @param xx numeric vector.
#' @param B number of bootstrap samples.
#' @param grid numeric vector used for evaluation
#' @return numeric vector with the density estimated over the given grid using BagHist
#'
#' @return estimation
#' @export
#' @import graphics
#' @examples
#' x <- rnorm(100)
#' BagHist(x, seq(-6, 6, length.out = 100), B = 50)
BagHist = function(xx, grid, B = 10) {
  n = length(xx)
  fin = 0
  for(i in 1:B){
    xb = xx[sample(n, replace = TRUE)]
    nbr = bropt(xb)$opt
    hs2 = hist(xb, breaks = mybreaks(xb, nbr), plot = FALSE, warn.unused = FALSE)
    fin = fin + predict_hist(hs2, grid)
  }
  fin / B
}
