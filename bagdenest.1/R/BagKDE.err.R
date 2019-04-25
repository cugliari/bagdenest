#' Title prevision and prevision error
#'
#' @param xx data vector
#' @param grille grid for density evaluation
#' @param B number of histograms to aggregate
#' @param dobs obervation
#'
#' @return  prevision and prevision error
#' @export
BagKDE.err = function(xx,grille=aa,B= 10,dobs) {
  # xx	data vector
  # grille	grid for density evaluation
  # B		number of histograms to aggregate
  # plerr	to plot error evolution
  fin = 0
  err00=NULL
  n = length(xx)
  for(i in 1:B) {
    xb = xx[sample(n,replace=T)]
    kk = kde(xb,h=bw.ucv(xb),eval.points=grille)
    fin = fin + kk$estimate
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #		if(i%%20 == 0) cat(i,">>")
  }
  #	cat("\n")
  list(prev=fin/B,erreur=err00)
}
