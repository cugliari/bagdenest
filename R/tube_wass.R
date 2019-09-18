#' @title  confidence band with histograms
#'
#' @description this function builds a confidence band from the histograms
#'
#' @param data ata of test sample and learning  (data$learn and data$test).
#'             test sample must be sorted  (sort(data$test)).
#' @param alpha parameter for qnorm
#' @param nbr number of breaks for histogram
#' @param plot  TRUE or FALSE, to plot confidence band
#' @param ... additional parameters for the histogram
#'
#' @return the lower edge and the upper edge of a confinement band
#'         and the aggregated version of the initial histogram.
#' @export
#'
tube_wass <- function(data, alpha = 0.1, nbr = 50, plot = FALSE, ...) {
  xb  <- data$train
  n   <- length(data$train)

  hs <- hist(xb, breaks = mybreaks(xb, nbr),
             warn.unused = FALSE, probability = TRUE, plot = plot, ...)

  m   <- length(hs$counts)
  c   <- qnorm(p = 1 - alpha / 2 / m) * sqrt(m / n)

  ln1    <- c(pmax(sqrt(hs$density) - c, 0)^2, 0)
  un1    <- c((sqrt(hs$density) + c)^2       , 0)

  hs2 <- predict_hist(hs, data$test)

  ln2    <- pmax(sqrt(hs2) - c, 0)^2
  un2    <- (sqrt(hs2) + c)^2

  if (plot) {
    lines(data$test, ln2, type = 's', col = 4, lwd = 2)
    lines(data$test, un2, type = 's', col = 4, lwd = 2)
  }

  return(list(punctual = hs2,
              tube = rbind(ln2, un2)))
}
