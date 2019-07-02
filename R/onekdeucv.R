#' @title kde estimation
#' @description this function estimates a density using kde
#' @param grid grid for density evaluation
#' @param xx data vector  for estimator bulding
#'
#' @return estimations
#' @export
onekdeucv = function(xx,grid) {
  ks::kde(xx,h=bw.ucv(xx),eval.points=grid)$estimate
}
