#' @title Construction of a variability band with frquency polygons.
#'
#' @description At each step, this function takes a new set of data, builds a
#'              frquency polygon,predicts and aggregates.
#' @param data data of test sample and learning  (data$learn and data$test).
#'             test sample must be sorted  (sort(data$test)).
#'
#' @param nbr number of sub-intervals of the initial partition
#' @param B Number of bootstrap replications
#' @param conf vector of the confidence threshold
#'
#' @return The vector of values of the aggregated version of a frquency polygon, a matrix of predictions
#'         of frquency polygon, matrix of the lower edge and the upper edge of a variability band.
#'
#' @export
tube_fp <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {

  xx <- data$train
  grille <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grille), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)

    mids <- hs2$mids
    delta <- diff(mybreaks)[1]

    newmids <- c(min(mids) - delta, mids, max(mids) + delta)

    mat[i, ] <- approxfun(x = newmids, y = c(0, hs2$density, 0))(grille)
  }
  mat[is.na(mat)] <- 0

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  return(list(punctual = punctual, bunch = mat, tube = tube))
}
