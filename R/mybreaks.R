#' Title breaks
#'
#' @param x data vector
#' @param nbr number of break sfor histogram
#'
#' @return  a breaks
#' @export
mybreaks = function(x = rnorm(100),nbr)
{
  mx = min(x)
  Mx = max(x)
  seq(mx,Mx, (Mx-mx)/nbr)
}
