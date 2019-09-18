#' @title Construction of a variability band with kernel density estimators.
#'
#' @description At each step, this function takes a new set of data, builds
#'              a kernel density estimator, predicts and aggregates.
#'
#' @param data ata of test sample and learning  (data$learn and data$test).
#'             test sample must be sorted  (sort(data$test)).
#' @param B Number of bootstrap replications
#' @param conf vector of the confidence threshold
#'
#' @return The vector of values of the aggregated version of a kernel density estimator, a matrix of predictions
#'         of kernel density estimator, matrix of the lower edge and the upper edge of a variability band.
#' @export
tube_KDE <- function(data, B = 10, conf = c(0.05, 0.95)) {

  xx <- data$train
  grid <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grid), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mat[i, ] <- onekdeucv(xx = xb, grid = grid)
  }

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  return(list(punctual = punctual, bunch = mat, tube = tube))
}


