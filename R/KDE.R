#' Title  prevision according to KDE
#'
#' @param dlearn  learning sample
#' @param dtest test sample
#' @param w ????
#' @param h step
#'
#' @return prvision
#' @export
KDE = function(dlearn,dtest,w=NULL,h=0.1){
  n = length(dlearn)
  if(is.null(w)) w=rep(1/n,n)
  yy = outer(dtest,dlearn,"-")/h
  aux = t(exp( -0.5 *(yy^2)) / sqrt(2* pi))
  colSums(w * aux) / h
}
