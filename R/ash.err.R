#' Title   prevision and prevision error
#'
#' @param xx data vector
#' @param aa grid for density evaluation
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#' @param dobs  observation
#'
#' @return  prevision and prevision error
#' @export
ash.err = function(xx,aa,nbr, B= 10,dobs) {
  mx=min(xx)
  Mx=max(xx)
  h=(Mx-mx)/nbr
  s=seq(mx,Mx,h)
  pred=0
  err0=NULL
  for (i in 1:B){
    hh=hist(xx,breaks=c(mx-0.5,s+(i-1)*h/B,Mx+0.5),plot=F,warn.unused = F)
    pred=pred+predict.hist(hh,sort(aa))
    predi=pred/i
    err0=rbind(err0,error(dobs,predi))
  }
  list(prev=predi,erreur=err0)
}
