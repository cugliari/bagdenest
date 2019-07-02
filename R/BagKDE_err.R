#' @title Density estimation with Bagkde_err
#'
#' @description  This function builds Kernel Densities estimators from the bootstrap samples.
#' Then averages  the estimates provided  by these estimators  and computes their estimation error.
#'
#' @param xx data vector  for estimator bulding
#' @param grid grid for density evaluation
#' @param dobs density values associated whit test sample
#' @param B number of histograms to aggregate
#'
#' @return  prevision and prevision error
#' @export
BagKDE_err = function(xx,grid,B= 10,dobs) {
  fin = 0
  err00=NULL
  n = length(xx)
  for(i in 1:B) {
    xb = xx[sample(n,replace=T)]
    kk =ks::kde(xb,h=bw.ucv(xb),eval.points=grid)
    fin = fin + kk$estimate
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
  }
  list(prev=fin/B,erreur=err00)
}
