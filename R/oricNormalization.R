#' REPLICATION ORIGIN BIAS normalization
#' This function normalizes the replication origin bias. Inspired by https://github.com/XiDsLab/CNV-BAC
#' 
#' @importFrom mgcv gam predict.gam
#' @param coverage dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param referenceLength
#' @param oriC_position position of oriC in the reference genomem can be one or multiple
#' @return A vector of normalized coverage values
#' @noRd
#' 
oricNormalization <- function(coverageVector, referenceLength, oriC_position){
  
  ## oriC position
  if(is.vector(oriC_position)){numOriC <- length(oriC_position)}else{numOriC <- length(oriC_position)}
  
  ## variables
  coverage_new <- rep(0, times=referenceLength) #prepare output vector
  RC_count <- data.frame(READEPTH=integer(),POS=integer(), oriDISTANCE=integer()) #reconstruct DF
  
  ################################################################################
  ## GET DISTANCE from ORIC
  RC_count[1:referenceLength,"POS"] <- seq(from = 1, to = referenceLength)
  RC_count[1:referenceLength,"READEPTH"] <- coverageVector
  
  # Circular genome correction and multiple oriC handling
  for(oriC in oriC_position){
    RC_count <- cbind(RC_count, abs(RC_count$POS-oriC))
    RC_count <- cbind(RC_count, abs(referenceLength - RC_count$POS + oriC))
    RC_count <- cbind(RC_count, referenceLength + RC_count$POS - oriC)
  }
  RC_count[,"oriDISTANCE"] <- apply(RC_count[,c(4:ncol(RC_count))],1,min)
  # positions in between two most distant oriC's are zeroed
  index <- intersect(which(RC_count$POS <= max(oriC)),which(RC_count$POS >= min(oriC)))
  RC_count[index,"oriDISTANCE"] <- 0
  
  #
  RC_count <- RC_count[,c(1:3)]
  
  # plot(RC_count$oriDISTANCE)
  
  ################################################################################

  ## BOOSTRAPING and PLOTTING
  XX <- RC_count[sample(nrow(RC_count), 50000, replace = FALSE), ]
  XX <- XX[order(XX$oriDISTANCE),]
  XX <- XX[XX$READEPTH<summary(abs(RC_count$READEPTH))[[5]],] #remove values behind 3rd quantile
  library(ggplot2)
  ggplot(aes(x=oriDISTANCE, y=READEPTH),data=XX) + geom_point()
  
  summary(RC_count[,"oriDISTANCE"])
  summary(RC_count[,"READEPTH"])
  
  ################################################################################
  
  ## CORRELATION BETWEEN DISTANCE and COVERAGE
  statistics <- cor.test(RC_count[,"oriDISTANCE"], RC_count[,"READEPTH"], method ="spearman", exact=FALSE)
  statistics$p.value
  statistics$estimate
  RHO <- statistics$estimate
  if (RHO >= 0.9){ #maybe different value
    # ## REPLACE OUTLIERS
    readdepthmedian <- median(abs(RC_count$READEPTH))
    readdepth_1Q <- summary(abs(RC_count$READEPTH))[[2]]
    readdepth_3Q <- summary(abs(RC_count$READEPTH))[[5]]
    readdepth_IQR <- IQR(abs(RC_count$READEPTH))
    #find outliers 1.5xIQR rule
    out_ind <- which(abs(RC_count$READEPTH)<(readdepth_3Q+1.5*readdepth_IQR) & abs(RC_count$READEPTH)>(readdepth_1Q-1.5*readdepth_IQR) )
    # replace outliers value with MEDIAN
    RC_count[out_ind,"READEPTH"] <- readdepthmedian
    # 
    # fit model with read-depth without outliers
    RC <- as.vector(RC_count$READEPTH)
    DIST <- as.vector(RC_count$oriDISTANCE)
    DIST <- scale(DIST)
    RC <- log2(RC)
    fit <- mgcv::bam(RC~s(DIST))
    RC <- as.data.frame(coverageVector)
    pred <- mgcv::predict.gam(fit,RC)
    expected <- RC_count$READEPTH + pred
    # RC_count$normalized <- expected
    coverageVector_new <- expected
    
    return(coverageVector_new)
  }
  else{
    return(coverageVector)
  }
}

