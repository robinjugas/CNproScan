#' Modified Z-score by Iglewicz and Hoaglin
#' This function loads a vector of coverage/read-depth values. It returns the indices of potential outliers. 
#' See https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
#' 
#' @param coverageVector Vector of coverage values, 3rd column of coverage DATAFRAME
#' @return A vector of outliers indices
#' @export
#' @noRd
#' 
mzscore <- function(coverageVector) {
  medianX <- median(coverageVector)
  MAD <- median(abs(coverageVector - medianX))
  M <- c()
  for (i in 1:length(coverageVector)) {
    M[i] <- 0.6745 * abs(coverageVector[i] - medianX) / MAD
  }
  idx <- which(M > 3.5)
  
  outliers <- coverageVector[idx]
  
  return(idx)
}