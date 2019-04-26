#' Title  prevision and prevision error
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#' @param dobs observation
#' @param alpha disturbance parameter
#'
#' @return prevision and prevision error
#' @export
rashfp.err = function(xx,grille=aa,nbr = 50, B= 10,dobs,alpha=1) {
  fin = 0
  err00=NULL
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    blabla=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(length(zz),0,alpha * blabla)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + approxfun(x=hs2$mids,y=hs2$density)(grille)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #cat("\n")
  list(prev=fin/B,erreur=err00[,1])
}
