#' GESD - Two-sided Generalized Extreme Studentized Deviate test
#' This function requires:
#' 1 - a dataframe with 2 columns - 1st column: position/coordinate 2nd column: coverage/read-depth 
#' 2- an upper bound for the suspected number of outliers (an integer)
#' It returns outliers indices as vector
#' See https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
#' 
#' @param coverageDF dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param estimatedOutliers Number of potential outliers from M-Zscore function
#' @return A vector of outliers indices from input
#' @export
#' @noRd
#' 
gesd <- function(coverageDF, estimatedOutliers) {
  # COMPUTATION
  coverageSignal <- coverageDF$COVERAGE
  alpha <- 0.05
  R <- replicate(estimatedOutliers, 0)
  Lambda <-  replicate(estimatedOutliers, 0)
  outliers <-  replicate(estimatedOutliers, 0)
  
  n <- length(coverageSignal)
  xCopy  <- coverageSignal
  outliersNumber <- 0
  
  for (i in 1:estimatedOutliers) {
    xMean <- mean(coverageSignal)
    xSD <- sd(coverageSignal)
    
    R[i] <- max(abs(coverageSignal - xMean) / xSD)
    
    # compute critical value
    p <- 1 - alpha / 2 / (n - i + 1)
    t <- qt(p, n - i - 1)
    Lambda[i] <-
      t * (n - i) / sqrt((n - i - 1 + t ^ 2) * (n - i + 1))
    
    if (R[i] > Lambda[i]) {
      outliersNumber <- i
    }
    
    idxMax <- which.max(abs(coverageSignal - xMean))
    outliers[i] <- coverageSignal[idxMax]
    coverageSignal <-  coverageSignal[-idxMax]
  }
  
  idx <- replicate(outliersNumber, 0)
  
  # match find values in original vector
  if (outliersNumber > 0){
    for (i in 1:outliersNumber) {
      idx[i] <- which(xCopy == outliers[i])[[1]]
      xCopy[idx[i]] <- NaN
      idx <- unlist(idx)
    }
  }
  
  # sort
  if (outliersNumber > 0) {
    coverageDF$isOutlier <- FALSE
    coverageDF$isOutlier[idx] <- TRUE
    return(as.vector(coverageDF$POS[idx]))
  } else {
    return()
    warning("No outliers found")
  }
}
