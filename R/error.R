#' @title Estimation error
#'
#' @description this function provides estimation  using quadratic error or Error log likelihood .
#'
#' @param obs  values associated with test sample.
#' @param prev   prevision
#'
#' @return quadratic error ,Error log likelihood
#' @export
error = function(obs,prev)	{
  er1 = mean((obs-prev)^2,na.rm=T)
  er2 = sum(ifelse(prev>0,log(prev),0))
  signif(c(er1,er2),3)
}
