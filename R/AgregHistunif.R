#' @title Aggregation of disturbed histograms.
#'
#' @description This function aggregates disturbed histograms.Each histogram is disturbed by adding
#'small random numbers to the sub-interval edges of the partitions from which they are constructed.
#'These small random numbers are generated according to a uniform law.
#'
#' @param xx data vector.this vector makes it possible to build histograms
#' @param grid  data vector.This vector makes it possible to provide estimates and estimation errors.
#' @param nbr number of breaks for histogram
#' @param B number of histograms to aggregate
#'
#' @return Estimation of the values of a density function .
#' @export
#' @import graphics
AgregHistunif = function(xx,grid,nbr = 50, B=10) {
  fin = 0
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
  }
  fin/B
}
