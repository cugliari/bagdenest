#' Title  prevision and prevision error
#'
#' @param xx  data vector
#' @param grille grid for density evaluation
#' @param B number of histograms to aggregate
#' @param dobs observation
#'
#' @return prevision and prevision error
#' @export
BagHistfp.err = function(xx,grille=aa, B= 10,dobs) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
  #
  n = length(xx)
  err00=NULL
  err002=NULL
  fin = 0
  fin2=0
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)    {
    xb = xx[sample(n,replace=T)]
    nbr=bropt(xb)$opt
    nbrfp=broptfp(xb)$opt
    hs=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    hs2=hist(xb,breaks=mybreaks(xb,nbrfp),plot=F,warn.unused = F)
    m <- hs2$mids
    h <- m[2] - m[1]
    m <- c(m[1] - h, m, m[length(m)] + h)
    d <- c(0, hs2$density, 0)
    fin  = fin  + approxfun(x = m, y = d, yright = 0, yleft = 0)(grille)

    fin2=fin2+ predict.hist(hs,grille)
    previ=fin/i
    previ2=fin2/i
    err00=rbind(err00,error(dobs,previ))
    err002=rbind(err002,error(dobs,previ2))
    #if(i%%20 == 0) cat(i,">>")
  }
  #  cat("\n")
  list(prev=fin/B,erreurbhfp=err00,erreurbh=err002)
}
