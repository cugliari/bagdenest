#' Title Prevision acording to Bagkde
#'
#' @param xx data vector
#' @param dobs  observation
#' @param B number of histograms to aggregate
#' @param grid  grid
#'
#' @return Prevision
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
