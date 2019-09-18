#' @title  confidence band with kernel density estimators
#' @description this function builds a confidence band from the kernel density estimators
#' @param data ata of test sample and learning  (data$learn and data$test).
#'             test sample must be sorted  (sort(data$test)).
#'
#' @return the lower edge and the upper edge of a confinement band
#'         and the aggregated version of the initial kernel density estimator.
#' @export

tube_KDEsm <- function(data) {# only at 90% (hard coded)
  fit <- sm::sm.density(data$train, h = bw.ucv(data$train), display = "none",
                        eval.grid = TRUE, eval.points = data$test)
  return(list(punctual = fit$estimate,
              tube = rbind(fit$lower, fit$upper)))
}
