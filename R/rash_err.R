#' @title density estimation and estimation error whith rash.
#' @description this function builds several histograms in disturbing by translation an initial partition
#'              of the histogram with a Gaussian noise and then aggregates them.
#' @param xx data vector
#' @param grid grid for density evaluation
#' @param nbr number of breaks for histogram
#' @param B number of histograms to aggregate
#' @param dobs density values associated whit test sample
#'
#' @return  estimation and estimattion error
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
