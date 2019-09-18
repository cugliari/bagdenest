#' @title Density estimation whit ash_err
#'
#' @description this function builds several histograms on partions obtained
#'              by non-randomly translating an initial partition.
#'
#' @param xx data vector for histograms construction .
#' @param aa grid for density evaluation
#' @param nbr number of break sfor histogram
#' @param B number of histograms to aggregate
#' @param dobs  density values associated whit test sample
#'
#' @return  prevision and prevision error
#' @export
ash_err = function(xx,aa,nbr, B= 10,dobs) {
  mx=min(xx)
  Mx=max(xx)
  h=(Mx-mx)/nbr
  s=seq(mx,Mx,h)
  pred=0
  err0=NULL
  for (i in 1:B){
    hh=hist(xx,breaks=c(mx-0.5,s+(i-1)*h/B,Mx+0.5),plot=F,warn.unused = F)
    pred=pred+predict_hist(hh,sort(aa))
    predi=pred/i
    err0=rbind(err0,error(dobs,predi))
  }
  list(prev=predi,erreur=err0)
}

