#' Title prevision and prevision error
#'
#' @param xx data vector
#' @param B number of histograms to aggregate
#' @param grid grid
#'
#' @return prevision and prevision error
#' @export
#' @import graphics
BagHist = function(xx,grid, B= 10) {
  n = length(xx)
  fin = 0
  for(i in 1:B)    {
    xb = xx[sample(n,replace=TRUE)]
    nbr=bropt(xb)$opt
    hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,grid)
  }
  fin/B
}
