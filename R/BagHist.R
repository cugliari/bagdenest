#' Title prevision and prevision error
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param B number of histograms to aggregate
#'
#' @return prevision and prevision error
#' @export
BagHist = function(xx,grille=aa, B= 10) {
  n = length(xx)
  fin = 0
  for(i in 1:B)    {
    xb = xx[sample(n,replace=TRUE)]
    nbr=bropt(xb)$opt
    hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    fin= fin + predicthist(hs2,grille)
  }
  fin/B
}
