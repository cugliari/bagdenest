#' Title prevision and prevision error
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param B  number of histograms to aggregate
#' @param dobs observation
#'
#' @return  prevision and prevision error
#' @export
BagHist.err = function(xx,grille=aa,B= 10,dobs) {
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
    fin= fin + predict.hist(hs2,grille)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #	cat("\n")
  list(prev=fin/B,erreur=err00)
}
