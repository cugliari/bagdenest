#' Title prevision and prevision error according to BagHist.err
#'
#' @param xx data vector
#' @param B  number of histograms to aggregate
#' @param dobs observation
#' @param grid grid
#'
#' @return  prevision and prevision error
#' @export
#' @import graphics
BagHist_err = function(xx,grid,B= 10,dobs) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
  #
  n = length(xx)
  err00=NULL
  fin = 0
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)    {
    xb = xx[sample(n,replace=T)]
    nbr=bropt(xb)$opt
    hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,grid)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #	cat("\n")
  list(prev=fin/B,erreur=err00)
}
