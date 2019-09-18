#'@title Density estimation with AgregHist.
#'
#'@description This function aggregates disturbed histograms.Each histogram is disturbed by adding
#'small random numbers to the sub-interval edges of the partitions from which they are constructed.
#'
#' @param xx data vector for histograms bulding.
#' @param nbr number of breaks  for histogram
#' @param B  number of histograms to aggregate
#' @param alpha  disturbance parametter
#' @param grid grid for density evaluation.
#'
#' @return Estimation.
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
