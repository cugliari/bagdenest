#' Title  kde estimation
#'
#' @param grid grid
#' @param xx data vector
#'
#' @return estimations
#' @export
onekdeucv = function(xx,grid) {
  # xx	data vector
  # grille	grid for density evaluation
  ks::kde(xx,h=bw.ucv(xx),eval.points=grid)$estimate
}
