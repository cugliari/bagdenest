#' @title Density estimation whit BagHist
#'
#' @description This function builds histograms from the bootstrap samples.
#' Then averages  the estimates provided  by these histograms and computes thes estimation errro.
#'
#' @param xx data vector for histograms construction .
#' @param B number of histograms to aggregate
#' @param grid grid for density evaluation
#' @return density estimated over the grid using bagged histogram
#'
#' @return estimation
#' @export
#' @import graphics
BagHist = function(xx,grid, B= 10) {
  n = length(xx)
  fin = 0
  for(i in 1:B){
    xb = xx[sample(n,replace=TRUE)]
    nbr=bropt(xb)$opt
    hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=FALSE,warn.unused = FALSE)
    fin= fin + predict_hist(hs2,grid)
  }
  fin/B
}
