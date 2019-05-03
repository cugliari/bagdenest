#' TitlePrevision acording to Bagkde
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param dobs  observation
#' @param B number of histograms to aggregate
#'
#' @return Prevision
#' @export
Bagkde <- function(xx, grille = aa, dobs, B = 10) {
  library(ks)
  n   = length(xx)
  kk=0
  for (i in 1:B){
    xb = xx[sample(n, replace = TRUE)]
    kk0 = kde(xb, h=bw.ucv(xb), eval.points = grille)$estimate
    kk=kk+kk0
  }
  kk/B
}
