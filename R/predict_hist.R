#' @title Prediction whith histogram
#' @description  this function
#'
#' @param hh  This function estimates values associated with a sample of data using the histogram.
#' @param x sample of observations
#'
#' @return prediction
#' @export
predict_hist = function(hh,x)	{
  res=NULL
  for(i in 1:length(x))
    res[i] = predict_hist_x(hh,x[i])
  res
}
