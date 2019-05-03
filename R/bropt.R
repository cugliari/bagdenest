#' Title needed for Hist & RASH
#'
#' @param x variable
#'
#' @return  a list
#' @export
bropt=function(x){
  Mgrid <- 2:(5 * floor(sqrt(length(x))))
  J     <- numeric(length(Mgrid))
  for(m in seq_along(Mgrid)) {
    J[m] <- riskhist(obs=x, Mgrid[m], xlim = c(min(x)-0.5, max(x)+0.5))
  }
  list(opt=Mgrid[which.min(J)])
}
