#' @title Density estimation whit BagHist_err
#'
#' @description This function builds histograms from the bootstrap samples.
#' Then averages  the estimates provided  by these histograms and computes thes estimation error.
#'
#' @param xx data vector for histograms construction .
#' @param B number of histograms to aggregate
#' @param grid grid for density evaluation
#' @param dobs density values associated whit test sample
#'
#' @return estimation and estimation error.
#' @export
#' @import graphics
BagHist_err = function(xx,grid,B= 10,dobs) {

  n = length(xx)
  err00=NULL
  fin = 0
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)    {
    xb = xx[sample(n,replace=T)]
    nbr=bropt(xb)$opt
    hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    fin= fin + predict_hist(hs2,grid)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
  }
  list(prev=fin/B,erreur=err00)
}
