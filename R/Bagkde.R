#' @title Density estimation with Bagkde
#'
#' @description  This function builds Kernel Densities estimators from the bootstrap samples.
#' Then averages  the estimates provided  by these estimators.
#'
#' @param xx data vector  for estimator bulding
#' @param grid grid for density evaluation
#' @param dobs density values associated whit test sample
#' @param B number of histograms to aggregate
#'
#' @return estimations
#' @export
Bagkde <- function(xx, grid, dobs, B = 10) {
  n   = length(xx)
  kk=0
  for (i in 1:B){
    xb = xx[sample(n, replace = TRUE)]
    kk0 =ks::kde(xb, h=bw.ucv(xb), eval.points = grid)$estimate
    kk=kk+kk0
  }
  kk/B
}
