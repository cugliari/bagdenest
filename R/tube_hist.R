#' @title Construction of a variability band using histograms.
#'
#' @description At each step, this function takes a new set of data, builds a
#'              histogram,predicts and aggregates.
#
#' @param data ata of test sample and learning  (data$learn and data$test).
#'             test sample must be sorted  (sort(data$test)).
#' @param nbr number of sub-intervals of the initial partition
#' @param B Number of bootstrap replication
#' @param conf vector of the confidence threshold
#'
#' @return  The vector of values of the aggregated version of a kernel density estimator, a matrix of predictions
#'          of kernel density estimator, matrix of the lower edge and the upper edge of a variability band.
#' @export
tube_hist <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {

  xx <- data$train
  grille <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grille), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)

    mat[i, ] <- predict_hist(hh = hs2, x = grille)
  }

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  return(list(punctual = punctual, bunch = mat, tube = tube))
}
