#' REPLICATION ORIGIN BIAS normalization
#' This function normalizes the replication origin bias. 
#' 
#' @importFrom mgcv gam predict.gam
#' @param coverage dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param oriC_position position of oriC in the reference genome
#' @return A vector of normalized coverage values
#' @noRd
#' 
oricNormalization <- function(coverage, oriC_position=1){
  ## 
  referenceLength <- length(coverage)
  coverage_new <- rep(0, times=referenceLength)
  RC_count <- data.frame(READEPTH=integer(),oriDISTANCE=integer())
  ## GET DISTANCE from ORIC
  RC_count[1:referenceLength,"READEPTH"] <- coverage
  RC_count[,"oriDISTANCE"] <- seq(from = 1, to = referenceLength)
  RC_count[,"oriDISTANCE"] <- abs(RC_count[,"oriDISTANCE"] - oriC_position) #distance from ORIC1
  # RC_count[,"oriDISTANCE"] <- referenceLength - RC_count[,"oriDISTANCE"] + oriC_position #distance from ORIC2
  # RC_count[,"oriDISTANCE"] <- referenceLength + RC_count[,"oriDISTANCE"] - oriC_position #distance from ORIC3
  
  ## BOOSTRAPING and PLOTTING
  # XX = RC_count[sample(nrow(RC_count), 50000, replace = FALSE), ]
  # library(ggplot2)
  # ggplot(aes(x=READEPTH, y=oriDISTANCE),data=XX) + geom_point()
  
  # summary(RC_count[,"oriDISTANCE"])
  # summary(RC_count[,"READEPTH"])
  
  ## CORRELATION BETWEEN DISTANCE and COVERAGE
  statistics <- cor.test(RC_count[,"oriDISTANCE"],RC_count[,"READEPTH"],method ="spearman",exact=FALSE)
  RHO <- statistics$estimate
  if (RHO >= 0.9){
    # ## REPLACE OUTLIERS
    # readdepthmedian <- median(abs(RC_count$READEPTH))
    # readdepth_1Q <- summary(abs(RC_count$READEPTH))[[2]]
    # readdepth_3Q <- summary(abs(RC_count$READEPTH))[[5]]
    # readdepth_IQR <- IQR(abs(RC_count$READEPTH))
    # #find outliers 1.5xIQR rule
    # out_ind <- which(abs(RC_count$READEPTH)<(readdepth_3Q+1.5*readdepth_IQR) & abs(RC_count$READEPTH)>(readdepth_1Q-1.5*readdepth_IQR) )
    # # replace outliers value with MEDIAN
    # RC_count[out_ind,"READEPTH"] <- readdepthmedian
    # 
    # fit model with read-depth without outliers
    RC <- RC_count$READEPTH
    DIST <- RC_count$oriDISTANCE
    fit <- mgcv::gam(RC~s(DIST))
    RC <- as.data.frame(coverage)
    pred <- mgcv::predict.gam(fit,RC)
    expected <- RC_count$READEPTH + pred
    # RC_count$normalized <- expected
    coverage_new <- expected
    
    return(coverage_new)
  }
  else{
    return(coverage)
  }
}

