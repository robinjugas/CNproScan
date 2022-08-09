#' cnvOutliers
#' This function detects outliers in whole coverage with GESD, merges close events together and returns a cell data type 
#' Inputs: coverage dataframe (1st col POS, 2nd col CCOVERAGE), path to fasta file
#' Output:  vector of modified COVERAGE values
#' 
#' @importFrom parallel splitIndices makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#' @param coverageDF dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param cores number of cores to be used by parallel package
#' @param peakDistanceThreshold distance between close peaks to be merged
#' @return list of vectors with CNV peaks coordinates 
#' @export
#' @noRd
cnvOutliers <- function(coverageDF,cores=2,peakDistanceThreshold=20){
  ## PARALLEL LOOP OUTLIERS DETECTION
  
  # cores <- parallel::detectCores(logical = FALSE) #should be odd?
  chunks <- parallel::splitIndices(length(coverageDF$COVERAGE), ncl = cores)
  
  #create the cluster
  my.cluster <- parallel::makeCluster(cores, type = "PSOCK")
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #  .packages = c("CNproScan"),
  Outliers <- foreach::foreach(i=1:cores, .export=c("gesd","mzscore"), .combine='c', .inorder=FALSE, .multicombine=FALSE) %dopar% {
    PARTIAL_DF <- coverageDF[chunks[[i]],]
    estimatedOutliers <- mzscore(PARTIAL_DF$COVERAGE)
    # if estimatedOutliers valid run GESD, check if is valid, else return empty vector
    if (!length(estimatedOutliers) == 0)
    {
      candidateOutliers <- gesd(PARTIAL_DF, length(estimatedOutliers))
      
      if (!length(candidateOutliers) == 0) {
        return(candidateOutliers)
      }
      else
      {
        candidateOutliers <- as.vector(c(0))
        return(candidateOutliers)
      }
    }
    else
    {
      candidateOutliers <- as.vector(c(0))
      return(candidateOutliers)
    }
    
  }
  
  #stop cluster
  parallel::stopCluster(cl = my.cluster)
  
  Outliers <- Outliers[Outliers!=0] #remove zeroes
  
  ##############################################################################
  ## MERGING CLOSE CNV EVENTS TOGETHER
  idxGESD <- sort(Outliers)
  differences <- diff(idxGESD)
  kk<-which(differences>1 & differences<peakDistanceThreshold)
  
  if(length(kk)){
    for (i in seq(1,length(kk))){
      CNVdistance <- idxGESD[kk[i]+1] - idxGESD[kk[i]]
      vec <- seq(idxGESD[kk[i]]+1, idxGESD[kk[i]]+CNVdistance-1) # coordinates to fill in the gap between 2 outliers
      #concatenate in between
      idxGESD <-c ( idxGESD[1:kk[i]] ,vec, idxGESD[(kk[i]+1):length(idxGESD)] )
      korekce <- CNVdistance-1 #correction - idxGESD longer by one
      kk <- kk+korekce # correction so for loop is not biased by increased kk
    }
  }
  ##############################################################################
  ## PREPROCESSING of CNV EVENTS into cell datatype
  idx <- idxGESD
  outliers_value <- coverageDF$COVERAGE[idx]
  differences <- diff(idx)
  k <- which(differences>1)
  
  if(length(k)){
    # https://stackoverflow.com/a/65758426 
    # list with a dimension attribute:
    peaks <- vector("list", length = length(k)+1)
    dim(peaks) <- c(length(k)+1, 1)
    for (i in seq(1,length(k)+1)){
      if (i==1){ #FIRST peak
        matrixtemp <- cbind(idx[1 : k[i]], outliers_value[1 : k[i]])
        peaks[[i, 1]] <- matrixtemp
      }else if(i==length(k)+1){ #LAST peak
        matrixtemp <- cbind(idx[ (k[i-1]+1) : length(idx)], outliers_value[ (k[i-1]+1) : length(outliers_value)])
        peaks[[i, 1]] <- matrixtemp
      }else{ #IN BETWEEN peaks
        matrixtemp <- cbind(idx[ (k[i-1]+1) : k[i] ], outliers_value[( k[i-1]+1) : k[i]  ])
        peaks[[i, 1]] <- matrixtemp
      }
    }
  }
  # if only single CNV is detected
  if(length(k)==0 & length(idxGESD)!=0)
  {
    peaks <- vector("list", length = 1)
    dim(peaks) <- c(1, 1)
    matrixtemp <- cbind(idx, outliers_value)
    peaks[[1, 1]] <- matrixtemp
  }
  
  return(peaks)
}

