#' @title  variability band
#' @description This function computes a variability band from a set of vectors
#' @param x data vector
#' @param punctual aggregation data vector
#' @param bunch estimations matrice
#' @param conf confidence threshold vector
#' @param plot plot validation
#'
#' @return the lower edge and the upper edge of a variability band
#' @export
baseboot <- function(x, punctual, bunch, conf = c(0.05, 0.95), plot = FALSE){

  err <- apply(bunch, 1, '-', punctual)
  CI  <- apply(err, 1, quantile, probs = conf)

  Linf <- punctual + CI[1, ]
  Lsup <- punctual + CI[2, ]

  if(plot) plotbunch(x, punctual, t(bunch), Linf, Lsup)

  return(rbind(Linf, Lsup))
}
