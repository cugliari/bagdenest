#'@title breaks
#'
#'@description This funtion  provides breaks for a density estimator
#'
#' @param x data vector
#' @param nbr number of breaks for an estimator
#'
#' @return  a breaks
#' @export
mybreaks = function(x = rnorm(100),nbr)
{
  mx = min(x)
  Mx = max(x)
  seq(mx,Mx, (Mx-mx)/nbr)
}
