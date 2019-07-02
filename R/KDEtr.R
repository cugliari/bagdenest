#' @title  estimation according to KDEtr
#'
#' @description This function makes density estimation whith Kernel Density Estimator using
#'              triangular kernel.
#'
#' @param dlearn  learning sample
#' @param dtest test sample
#' @param w  a vector of 1/n , where n is the length learning sample .
#' @param h step
#'
#' @return Estimation
#' @export
KDEtr = function(dlearn,dtest,w=NULL,h=0.1)
{
  n = length(dlearn)
  if(is.null(w)) w=rep(1/n,n)
  yy = outer(dtest,dlearn,"-")/h
  noyau = t((1- abs(yy))*ind2(yy,-1,1) )
  colSums(w * noyau) / h
}
