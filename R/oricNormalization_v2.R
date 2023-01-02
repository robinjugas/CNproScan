#' REPLICATION ORIGIN BIAS normalization  v2
#' This function normalizes the replication origin bias. Inspired by Yoon GC normalization
#' 
#' @param coverageVector just vector of coverage 1:referenceLength
#' @param referenceLength length of reference chromosome
#' @param oriC_position position of oriC in the reference genomem can be one or multiple
#' @return A vector of normalized coverage values
#' @export
#' @noRd
oricNormalization_v2 <- function(coverageVector, referenceLength, oriC_position){
  ##
  coverageVECTOR_new <- coverageVector
  windowLength <- 100
  level_round <- -4 #-5 10^5 / -4 10^4
  ## oriC position input
  if(is.vector(oriC_position)){numOriC <- length(oriC_position)}else{numOriC <- length(oriC_position)}
  
  ################################################################################
  ## remove outliers first for higher robustness
  # 1.5 IQR rule
  coverage_1Q <- summary(coverageVector)[[2]]
  coverage_3Q <- summary(coverageVector)[[5]]
  coverage_IQR <- IQR(coverageVector)
  
  if((coverage_1Q-1.5*coverage_IQR) > min(coverageVector)){
    extremesLow <- which(coverageVector<=(coverage_1Q-1.5*coverage_IQR))
  }else(extremesLow <- c())
  
  if((coverage_3Q+1.5*coverage_IQR) < max(coverageVector)){
    extremesHigh <- which(coverageVector>=(coverage_3Q+1.5*coverage_IQR))
  }else(extremesHigh <- c())
  
  coverageVector[extremesLow] <- median(coverageVector)
  coverageVector[extremesHigh] <- median(coverageVector)
  
  ################################################################################
  ## DATA.FRAME with windows
  WINDOWS <- data.frame(START=integer(), STOP=integer(),POS=integer(),READDEPTH=integer(),oriDISTANCE=integer())
  WINDOWS[1:length(seq(from=1, to=referenceLength, by=windowLength)),"START"] <- seq(from=1, to=referenceLength, by=windowLength)
  WINDOWS$STOP <- WINDOWS$START+windowLength-1
  WINDOWS <- WINDOWS[-nrow(WINDOWS)] #delete last window
  WINDOWS$STOP[nrow(WINDOWS)] <- referenceLength #extend the last windows to the end
  WINDOWS$POS <- WINDOWS$START+round((WINDOWS$STOP-WINDOWS$START)/2)
  
  ################################################################################
  ## Calculate oriC distance and READ count in window
  RC_count <- rep(0, times=nrow(WINDOWS))
  for (i in 1:nrow(WINDOWS)){
    RC_count[i] <- mean(coverageVector[WINDOWS[i,"START"]:WINDOWS[i,"STOP"]]) # average read depth in windows
  }
  WINDOWS[, "READDEPTH"] <- round(RC_count)
  ## median read depth of windows
  medReadCount <- median(RC_count) #median RC for all bins
  
  ################################################################################
  ## Circular genome correction and multiple oriC handling
  ORICDISTANCE <- rep(0, times=nrow(WINDOWS))
  for(oriC in oriC_position){
    # distance to oriC
    ORICDISTANCE <- cbind(ORICDISTANCE, abs(WINDOWS$POS-oriC))
    ORICDISTANCE <- cbind(ORICDISTANCE, abs(referenceLength - WINDOWS$POS + oriC))
    ORICDISTANCE <- cbind(ORICDISTANCE, referenceLength + WINDOWS$POS - oriC)
  }
  
  WINDOWS[,"oriDISTANCE"] <- apply(ORICDISTANCE[,c(2:ncol(ORICDISTANCE))],1,min)
  WINDOWS[,"oriDISTANCE"] <- round(WINDOWS[,"oriDISTANCE"],level_round)
  
  # positions in between two most distant oriC's are zeroed ? if oriC are defined near beginning and end it makes whole genome zero oriCdistance
  # index <- intersect(which(WINDOWS$POS <= max(oriC_position)),which(WINDOWS$POS >= min(oriC_position)))
  # WINDOWS[index,"oriDISTANCE"] <- 0
  
  ################################################################################
  ## oriC dataframe oriC_tab - 1.col - distance to ORIC in thousands 2.col - median coverageVECTOR of all windows with given distance to ORIC
  oriC_tab <- data.frame(oriDISTANCE=integer(), MED_COV=integer())
  oriC_tab[1:length(unique(WINDOWS$oriDISTANCE)), "oriDISTANCE"] <- unique(WINDOWS$oriDISTANCE)
  oriC_tab[, "MED_COV"] <- 0
  
  # get MED COV for given distance to ORIC
  for(i in 1:nrow(oriC_tab)){
    position <- which(WINDOWS$oriDISTANCE==oriC_tab[i, "oriDISTANCE"])
    oriC_tab[i, "MED_COV"]  <- median(na.omit(WINDOWS$READDEPTH[position]))
  }
  
  ## find corresponding median read-depth of given GC percentage score for every region
  for  (i in 1:nrow(WINDOWS)){
    idx <- which(oriC_tab[, "oriDISTANCE"] == WINDOWS[i, "oriDISTANCE"])
    WINDOWS[i, "MEDIANRC"] <- oriC_tab[idx, "MED_COV"]
  }
  
  
  ################################################################################
  ## CORRELATION BETWEEN DISTANCE and COVERAGE
  #  temp XX
  XX <- as.data.table(oriC_tab)
  XX <- unique(XX,by="oriDISTANCE")
  XX <- unique(XX,by="MED_COV")
  XX <- XX[order(XX$oriDISTANCE),]

  statistics <- cor.test(XX$oriDISTANCE, XX$MED_COV, method ="spearman", exact=TRUE)
  # statistics$p.value<0.05
  # statistics$p.value
  # statistics$estimate
  
  RHO <- statistics$p.value
  
  if (RHO <= 0.05){ #maybe different value?
   
    
    ## oriC Normalization
    for  (i in 1:nrow(WINDOWS)){
      start <- WINDOWS[i, "START"]
      stop <- WINDOWS[i, "STOP"]
      
      if(WINDOWS[i, "MEDIANRC"]>1){
        coeficient <- medReadCount/WINDOWS[i, "MEDIANRC"]
        coverageVECTOR_new[start:stop] <- coverageVECTOR_new[start:stop]*coeficient
      }# no normalization for extreme values
    }
    return(coverageVECTOR_new)
    
  } else{
    return(coverageVECTOR_new)
  }
  
}

