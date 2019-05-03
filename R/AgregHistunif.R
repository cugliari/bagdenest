#' Title  ???????
#'
#' @param xx  data vector
#' @param grille grid for density evaluation
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#'
#' @return ??????
#' @export
AgregHistunif = function(xx,grille=aa,nbr = 50, B=10) {
  # xx	data vector
  # grille	grid for density evaluation
  # br	number of break sfor histogram
  # B		number of histograms to aggregate
  #Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
  fin = 0
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
  z=diff(zz)
  h=z[1]
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    #blabla=sqrt(abs(min(diff(zz))))
    newb = zz + runif(length(zz),0,h)
    #newb=sort(newb)
    #if(min(newb) > mx) newb= c(mx,newb)
    #if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=c(mx,newb,Mx),plot=F,warn.unused = F)
    fin= fin + predicthist(hs2,grille)
    #if(i%%20 == 0) cat(i,">>")
  }
  #cat("\n")
  fin/B
}
