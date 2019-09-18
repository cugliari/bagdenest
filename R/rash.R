#' @title Density estimation rash.
#' @description This function builds several histograms in disturbing by translation an intial partition of the histogram
#'with a Gaussian noise and then aggregates them.Each partition on which a histogram is built is a
#'translation of the initial partition by adding this Gaussian noise to the edges of its sub-intervals.
#'
#' @param xx  data vector
#' @param nbr number of breaks for histogram
#' @param B number of histograms to aggregate
#' @param grid grid for density evaluation
#'
#' @return Prevision
#' @export
rash = function(xx,grid,nbr = 50, B=10) {
  fin = 0
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
  }
  fin/B
}
