#' @title  estimation according to KDE
#'
#' @description This function makes density estimation whith Kernel Density Estimator using
#'              gaussian kernel.
#'
#' @param dlearn  learning sample
#' @param dtest test sample
#' @param w  a vector of 1/n , where n is the length learning sample .
#' @param h step
#'
#' @return Estimation
#' @export
KDE = function(dlearn,dtest,w=NULL,h=0.1){
  n = length(dlearn)
  if(is.null(w)) w=rep(1/n,n)
  yy = outer(dtest,dlearn,"-")/h
  aux = t(exp( -0.5 *(yy^2)) / sqrt(2* pi))
  colSums(w * aux) / h
}
