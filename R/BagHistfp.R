#' @title Density Estimation with Bagged Frequency Polygons (BagHistFP)
#'
#' @description Estimate the density of a univariate data vector using BagHistFP. First,
#' obtain bootstrap samples of the input vector data. For each sample, a histogram is
#' computed using an optimized number of break points. On top of it, a frequency polygon
#' is obtained. The BagHistFP estimation is the mean of the individual frequency polygons.
#'
#' @param xx numeric vector.
#' @param B number of bootstrap samples.
#' @param grid numeric vector used for evaluation
#' @return numeric vector with the density estimated over the given grid using BagHistFP
#'
#' @export
#' @import graphics
#' @examples
#' x <- rnorm(100)
#' BagHistfp(x, seq(-6, 6, length.out = 100), B = 50)

BagHistfp = function(xx, grid, B = 10) {
  n = length(xx)
  fin = 0
  fin2 = 0
  for (i in 1:B)    {
    xb = xx[sample(n, replace = TRUE)]
    nbr = bropt(xb)$opt
    nbrfp = broptfp(xb)$opt
    hs = hist(
      xb,
      breaks = mybreaks(xb, nbr),
      plot = FALSE,
      warn.unused = FALSE
    )
    hs2 = hist(
      xb,
      breaks = mybreaks(xb, nbrfp),
      plot = FALSE,
      warn.unused = FALSE
    )
    m <- hs2$mids
    h <- m[2] - m[1]
    m <- c(m[1] - h, m, m[length(m)] + h)
    d <- c(0, hs2$density, 0)
    fin  = fin  + approxfun(
      x = m,
      y = d,
      yright = 0,
      yleft = 0
    )(grid)
    fin2 = fin2 + predict_hist(hs, grid)
  }

  list(bhfp = fin / B, bh = fin2 / B)
}
