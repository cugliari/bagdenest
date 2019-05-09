#' Title Prediction according to the histogram
#'
#' @param hh ????
#' @param x ?????
#'
#' @return prediction
#' @export
predict_hist = function(hh,x)	{
  res=NULL
  for(i in 1:length(x))
    res[i] = predict_hist_x(hh,x[i])
  res
}
