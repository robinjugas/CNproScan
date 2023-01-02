#' REPLICATION ORIGIN BIAS normalization
#' This function normalizes the replication origin bias. Credits to: https://github.com/XiDsLab/CNV-BAC
#' 
#' @importFrom mgcv gam predict.gam
#' @param coverage dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param referenceLength length of reference chromosome
#' @param oriC_position position of oriC in the reference genomem can be one or multiple
#' @return A vector of normalized coverage values
#' @noRd
#' 
oricNormalization_GAM <- function(coverageVector, referenceLength, oriC_position){
  
  ## oriC position
  if(is.vector(oriC_position)){numOriC <- length(oriC_position)}else{numOriC <- length(oriC_position)}
  
  ## variables
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
  
  ################################################################################
  
  # # plot(RC_count$oriDISTANCE)
  # summary(coverageVector)
  # summary(RC_count$READEPTH)
  # 
  # ## BOOSTRAPING and PLOTTING
  # XX <- RC_count[sample(nrow(RC_count), 50000, replace = FALSE), ]
  # XX <- XX[order(XX$oriDISTANCE),]
  # XX <- XX[XX$READEPTH<summary(abs(RC_count$READEPTH))[[5]],] #remove values behind 3rd quantile
  # library(ggplot2)
  # ggplot(aes(x=oriDISTANCE, y=READEPTH),data=XX) + geom_point()
  # 
  # summary(RC_count[,"oriDISTANCE"])
  # summary(RC_count[,"READEPTH"])
  # summary(RC)
  # summary(oriDIST)
  ################################################################################
  
  ## CORRELATION BETWEEN DISTANCE and COVERAGE
  statistics <- cor.test(RC_count[,"oriDISTANCE"], RC_count[,"READEPTH"], method ="spearman", exact=FALSE)
  statistics$p.value
  statistics$estimate
  RHO <- statistics$estimate
  
  if (RHO >= 0.5){ #maybe different value?
  
    # remove outlier residuals
    RC <- as.vector(RC_count$READEPTH)
    oriDIST <- as.vector(RC_count$oriDISTANCE)
    
    fit <- mgcv::gam(RC~s(oriDIST))
    
    residuals <- fit$residuals
    iqr <- IQR(residuals)
    index <- which(residuals > quantile(residuals,0.25)-1.5*iqr & residuals < quantile(residuals,0.75)+1.5*iqr)
    
    # fit model I
    RCnew <- RC[index]
    oriDISTnew <- oriDIST[index]
    
    fit <- mgcv::gam(RCnew~s(oriDISTnew))
    
    DF <- data.frame(oriDISTnew=oriDIST)
    pred <- mgcv::predict.gam(fit,DF)
    
    RCexpected <- as.vector(RC + pred)
    IQR <- quantile(RCexpected,0.75)-quantile(RCexpected,0.25)
    index2 <- which(RCexpected > quantile(RCexpected,0.25)-1.5*IQR & RCexpected < quantile(RCexpected,0.75)+1.5*IQR)
    
    # fit model II
    RCnew <- RC[index2]
    oriDISTnew <- oriDIST[index2]
    
    fit <- mgcv::gam(RCnew~s(oriDISTnew))
    DF <- data.frame(oriDISTnew=oriDIST)
    pred <- mgcv::predict.gam(fit,DF)
    
    RCexpected <- as.vector(coverageVector + pred)
    coverageVector_new <- as.integer(RCexpected)
    return(coverageVector_new)
    
  }
  else{
    return(coverageVector)
  }
}

# ##COMPARISON
# summary(coverageVector)
# summary(coverageVector_new)
# diff <- abs(coverageVector_new-coverageVector)
# 
# # plot genomes
# par( mfrow= c(3,1))
# plot(coverageVector)
# plot(coverageVector_new)
# plot(diff)