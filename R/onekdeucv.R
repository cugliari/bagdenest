#' Title ??
#'
#' @param xx data vector
#' @param grille grille	grid for density evaluation
#'
#' @return ??
#' @export
onekdeucv = function(xx,grille=aa) {
  # xx	data vector
  # grille	grid for density evaluation
  kde(xx,h=bw.ucv(xx),eval.points=grille)$estimate
}
