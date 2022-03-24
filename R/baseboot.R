#' @title Variability band
#' @description Computation of a variability band using base bootstrap
#' @param x data vector
#' @param pointwise a vector with the pointwise density estimation
#' @param bunch bootstrap replicated of the estimation density (matrix)
#' @param conf confidence quantiles
#' @param plot logical, make a validation plot
#'
#' @return matrix with lower and upper bounds of the variability band
#' @export
baseboot <- function(x, pointwise, bunch, conf = c(0.05, 0.95), plot = FALSE){

  err <- apply(bunch, 1, '-', pointwise)
  CI  <- apply(err, 1, quantile, probs = conf)

  Linf <- pointwise + CI[1, ]
  Lsup <- pointwise + CI[2, ]

  if(plot) plotbunch(x, pointwise, t(bunch), Linf, Lsup)

  return(rbind(Linf, Lsup))
}
