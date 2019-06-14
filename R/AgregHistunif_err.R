#'@title Aggregation of disturbed histograms and estimation error.
#'
#'@description This function aggregated histograms disturbed.each histogram is disturbed with a uniform variable.
#' Thus for each histogram a simulation of this uniform variable is performed.
#'
#' @param xx data vector.this vector makes it possible to build histograms
#' @param grid  data vector.This vector makes it possible to provide estimates and estimation errors.
#' @param nbr number of breaks for histogram
#' @param B number of histograms to aggregate
#' @param dobs density values associated with test sample.
#'
#' @return Estimation of the values of a density function And estimation error.
#' @export
AgregHistunif_err = function(xx,grid,nbr = 50, B= 10,dobs) {
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
    fin= fin + predict_hist(hs2,grid)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #cat("\n")
  list(prev=fin/B,erreur=err00)
}
