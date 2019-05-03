#' Title  ???????
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param nbr
#' @param B  number of histograms to aggregate
#' @param alpha  disturbance parametter
#'
#' @return ????
#' @export
AgregHist = function(xx,grille=aa,nbr = 50, B=10, alpha=1) {
  fin = 0
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    blabla=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(length(zz),0,alpha*blabla)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + predicthist(hs2,grille)
  }
  fin/B
}
