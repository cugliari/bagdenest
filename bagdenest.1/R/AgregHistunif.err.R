#' Title prevision and prevision error
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#' @param dobs observation
#'
#' @return
#' @export
#'
#' @examples
AgregHistunif.err = function(xx,grille=aa,nbr = 50, B= 10,dobs) {
  fin = 0
  err00=NULL
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  z=diff(zz)
  h=z[1]
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    newb = zz + runif(length(zz),0,h)
    hs2=hist(xx,breaks=c(mx,newb,Mx),plot=F,warn.unused = F)
    fin= fin + predict.hist(hs2,grille)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #cat("\n")
  list(prev=fin/B,erreur=err00)
}
