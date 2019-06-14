#'@title Aggregation of disturbed histograms.
#'
#'@description This function aggregates disturbed histograms.Each histogram is disturbed by adding to each of its estimates an observation of a centered normal variable
#' whose variance is coeficiant "alpha" which is the pertubation parameter.For each histogram a simulation of this normal variable is performed.
#'
#' @param xx data vector for histograms bulding.
#' @param nbr number of breaks  for histogram
#' @param B  number of histograms to aggregate
#' @param alpha  disturbance parametter
#' @param grid  data vector.this vector makes it possible to provide estimations and estimation errors.
#'
#' @return Estimation of the values of a density function.
#' @export
#' @import graphics
AgregHist = function(xx,grid,nbr = 50, B=10, alpha=1) {
  fin = 0
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    eps=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(length(zz),0,alpha*eps)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,grid)
  }
  fin/B
}
