#' @title Aggregation of disturbed histograms and estimation error.
#'
#' @description This function aggregates disturbed histograms an provides estimation error from the aggregation.
#' This function aggregates disturbed histograms.Each histogram is disturbed by adding small random numbers to
#' the sub-interval edges of the partitions from which they are constructed.These small random numbers are
#' generated according to a normal law.
#'
#' @param xx data vector . this vector makes it possible to build histograms.
#' @param grid  data vector.this vector makes it possible to provide estimates and estimation errors.
#' @param nbr number of breaks for histogram
#' @param B number of histograms to aggregate
#' @param dobs density values associated with test sample.
#' @param alpha pertubation parametter
#'
#' @return  Estimation of the values of a density function And estimation error.
#' @export
#' @import graphics
AgregHist_err = function(xx,grid,nbr = 50, B= 10,dobs,alpha=1) {
  fin = 0
  err00=NULL
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    eps=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(length(zz),0,alpha * eps)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,eps)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
  }
  list(prev=fin/B,erreur=err00)
}

