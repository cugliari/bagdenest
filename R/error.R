#' Title ????
#'
#' @param obs  observation
#' @param prev   Pprevision
#'
#' @return ?????
#' @export
error = function(obs,prev)	{
  #1) Erreur quadratique (c'est le ISE).
  er1 = mean((obs-prev)^2,na.rm=T)
  #2) Erreur log vraissemblance
  er2 = sum(ifelse(prev>0,log(prev),0))
  #2bis) Distance KL
  #	er22 = mean(ifelse(obs>0,log(obs),0)) - mean(ifelse(prev>0,log(prev),0))
  #3) Erreur de dissimilarit? d'Hellinger
  #	er3 =mean((sqrt(obs)-sqrt(prev))^2)
  signif(c(er1,er2),3)
}
