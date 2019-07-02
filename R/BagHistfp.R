#' @title Density estimation whit BagHistfp
#'
#' @description This function builds "Frequency Folygones"  from the bootstrap samples.
#' Then averages  the estimates provided  by these estimators .
#'
#' @param xx data vector for Frequency Folygone construction .
#' @param B number of Frequency Folygone to aggregate
#' @param grid grid for density evaluation
#'
#' @return estimation
#' @export
#' @import graphics
BagHistfp = function(xx,grid, B= 10) {

  n = length(xx)
  fin = 0
  fin2=0
  for(i in 1:B)    {
    xb = xx[sample(n,replace=TRUE)]
    nbr=bropt(xb)$opt
    nbrfp=broptfp(xb)$opt
    hs=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    hs2=hist(xb,breaks=mybreaks(xb,nbrfp),plot=F,warn.unused = F)
    m <- hs2$mids
    h <- m[2] - m[1]
    m <- c(m[1] - h, m, m[length(m)] + h)
    d <- c(0, hs2$density, 0)
    fin  = fin  + approxfun(x = m, y = d, yright = 0, yleft = 0)(grid)
    fin2 = fin2 + predict_hist(hs,grid)
  }

  list(bhfp=fin/B,bh=fin2/B)
}
