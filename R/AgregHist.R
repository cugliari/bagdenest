#' Title  prevision and prevision error according to AgregHist
#'
#' @param xx data vector
#' @param nbr number of break sfor histogram
#' @param B  number of histograms to aggregate
#' @param alpha  disturbance parametter
#' @param grid  grid
#'
#' @return prevision
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
