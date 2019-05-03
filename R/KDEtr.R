#' Title ???????
#'
#' @param dlearn  learning saple
#' @param dtest test sample
#' @param w ?????
#' @param h ?????
#'
#' @return ?????
#' @export
KDEtr = function(dlearn,dtest,w=NULL,h=0.1)
{
  library(ks)
  n = length(dlearn)
  if(is.null(w)) w=rep(1/n,n)
  yy = outer(dtest,dlearn,"-")/h
  noyau = t((1- abs(yy))*ind2(yy,-1,1) )
  colSums(w * noyau) / h
}
