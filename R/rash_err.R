#' Title  prevision and prevision error
#'
#' @param xx data vector
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#' @param dobs observation
#' @param grid  grid
#'
#' @return  prevision and prevision error
#' @export
#' @import graphics
rash_err = function(xx,grid,nbr = 50, B= 10,dobs) {
  fin = 0
  err00=NULL
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    eps=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(1,0,eps)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,grid)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #		if(i%%20 == 0) cat(i,">>")
  }
  #	cat("\n")
  list(prev=fin/B,erreur=err00)
}
