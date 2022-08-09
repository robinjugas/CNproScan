#' cnvDeletions
#' This function detects zero coverage and send that into a dataframe of deletions. Also merging close deletions together. 
#' 
#' @param coverageDF dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param peakDistanceThreshold distance between close peaks to be merged
#' @return Dataframe of detected deletions DEL_DF
#' @export
#' @noRd
#' 
cnvDeletions <- function(coverageDF, peakDistanceThreshold=20){
  ################################################################################
  ## DELETING crossovers at the start and at the end and replace them by mean of whole COVERAGE
  cutoff <- 100
  coverageDF$COVERAGE[1:cutoff] <- rep(mean(coverageDF$COVERAGE), times=length(coverageDF$COVERAGE[1:cutoff]))
  end <- length(coverageDF$COVERAGE)
  coverageDF$COVERAGE[(end-cutoff):end] <- rep(mean(coverageDF$COVERAGE), times=length(coverageDF$COVERAGE[(end-cutoff):end]))
  
  ################################################################################
  ## ZERO COVERAGE -> DELETIONS 
  deletions <- which(coverageDF$COVERAGE == 0)
  
  ################################################################################
  ## MERGING CLOSE DELETIONS TOGETHER
  idxDEL <- sort(deletions)
  differences <- diff(idxDEL)
  kk <- which(differences > 1 & differences < peakDistanceThreshold)
  
  if(length(kk)){
    for (i in seq(1,length(kk))){
      CNVdistance <- idxDEL[kk[i]+1] - idxDEL[kk[i]]
      vec <- seq(idxDEL[kk[i]]+1, idxDEL[kk[i]]+CNVdistance-1) # coordinates to fill in the gap between 2 outliers
      idxDEL <-c ( idxDEL[1:kk[i]] ,vec, idxDEL[(kk[i]+1):length(idxDEL)] ) #concatenate in between
      korekce <- CNVdistance-1 #correction - idx longer by one
      kk <- kk+korekce # correction so for loop is not biased by increased kk
    }
  }
  
  ################################################################################
  ## PREPROCESSING of DELETIONS into cell datatype
  idx <- idxDEL
  outliers_value <- coverageDF$COVERAGE[idx]
  differences <- diff(idx)
  k <- which(differences>1)
  
  if(length(k)){
    ## https://stackoverflow.com/a/65758426 
    ## list with a dimension attribute:
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
    ## SAVE into dataframe
    DEL_DF <- data.frame(ID=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), TYPE=character())
    for(i in 1:length(peaks)){
      start <-  peaks[[i, 1]] [1, 1]
      nrows <- nrow(peaks[[i, 1]])
      stop <- peaks[[i, 1]] [nrows, 1]
      DEL_DF[i,"ID"] <- i
      DEL_DF[i,"START"] <- start
      DEL_DF[i,"END"] <- stop
      DEL_DF[i,"LENGTH"] <- (stop - start)
      DEL_DF[i,"COVERAGE"] <- mean(coverageDF$COVERAGE[start:stop])
      DEL_DF[i,"TYPE"] <- "DEL"
    }
    ## DELETIONS COVERAGE VALUE REPLACED BY MEAN - INFLUENCES GESD OUTLIERS DETECTION
    coverageDF$COVERAGE[idxDEL] <- mean(coverageDF$COVERAGE)
  }
  ## if only single deletion is detected
  if(length(k) == 0 & length(idxDEL) != 0){
    DEL_DF <- data.frame(ID=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), TYPE=character())
    start <-  idxDEL[1]
    stop <- idxDEL[length(idxDEL)]
    DEL_DF[1, "ID"] <- 1
    DEL_DF[1, "START"] <- start
    DEL_DF[1, "END"] <- stop
    DEL_DF[1, "LENGTH"] <- (stop - start)
    DEL_DF[1, "COVERAGE"] <- mean(coverageDF$COVERAGE[start:stop])
    DEL_DF[1, "TYPE"] <- "DEL"
  }
  
  return(list(coverageDF, DEL_DF))
}

