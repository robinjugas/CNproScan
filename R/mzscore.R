#' Modified Z-score by Iglewicz and Hoaglin
#' This function loads a vector of coverage/read-depth values. It returns the indices of potential outliers. 
#' See https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
#' 
#' @param infile Path to the input file
#' @return A vector of outliers indices from input
#' @export
mzscore <- function(x) {
  medianX = median(x)
  MAD = median(abs(x - medianX))
  M <- c()
  for (i in 1:length(x)) {
    M[i] <- 0.6745 * abs(x[i] - medianX) / MAD
  }
  idx <- which(M > 3.5)
  outliers <- x[idx]
  return(outliers)
}