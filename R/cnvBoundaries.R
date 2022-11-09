#' cnvBoundaries
#' Formats outliers into CNV DUPLICATIONS dataframe and detects more precise CNV boundaries based on slope of a peak. 
#' 
#' @param peaksCELL list of vectors with CNV peaks coordinates 
#' @param coverageDF dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param averageCoverage value of average coverage
#' @param step number of bases to calculate slope
#' @param TRIGGER number of changes of slope derivation
#' @return CNV dataframe
#' @export
#' @noRd
#' 
cnvBoundaries <- function(peaksCELL, coverageDF, averageCoverage, step=11, TRIGGER=5){
  ## CNV EVENTS BORDERS DETECTION based on slope of a peak
  ## works only on POSITIVE PEAKS, NOT NEGATIVE PEAKS - DELETIONS ?
  peaksPolished <- vector("list", length=length(peaksCELL))
  dim(peaksPolished) <- c(length(peaksCELL), 1)
  
  CNVcoverage <- c()
  CNVcoverageRelative <- c()
  
  coverageSignal <- coverageDF$COVERAGE
  
  if(length(peaksCELL)>0){
    for (i in seq(1, length(peaksCELL))) {
      start <-  peaksCELL[[i, 1]][1, 1]
      nrows <- nrow(peaksCELL[[i, 1]])
      stop <- peaksCELL[[i, 1]] [nrows, 1]
      
      avg_cov_peak <- mean(coverageSignal[start:stop]) #CNVcoverage
      
      
      if (avg_cov_peak <= averageCoverage) {
        #DELETIONS \/ -> SKIP
        #SLOPE \ peak start
        count1 <- 0
        trigger1 <- 0
        # while (trigger1 < TRIGGER) {
        #   if (start <= 1 | (start - count1 * step - step) <= 1)
        #   {
        #     break
        #   }
        #   else
        #   {
        #     vektor <- coverageSignal[(start - count1 * step - step):(start - count1 * step)]
        #     # mod_vektor <- mode(vektor) #matlab
        #     mod_vektor <- sort(table(vektor), decreasing = FALSE)[1]
        #   }
        #   if (any(vektor == 0)) {
        #     break
        #   } #if zero in vector, ends - cause zero coverage
        #   
        #   count1 <- count1 + 1
        #   if (mod_vektor >= avg_cov_peak)
        #   {
        #     trigger1 <- trigger1 + 1
        #   }
        # }
        # SLOPE / peak stop
        count2 <- 0
        trigger2 <- 0
        # while (trigger2 < TRIGGER) {
        #   if (stop >= length(coverageSignal) |
        #       (stop + count2 * step + step) >= length(coverageSignal))
        #   {
        #     break
        #   }
        #   else{
        #     vektor <- coverageSignal[(stop + count2 * step):(stop + count2 * step + step)]
        #     
        #     # mod_vektor <- mode(vektor); #matlab
        #     mod_vektor <- sort(table(vektor), decreasing = FALSE)[1]
        #   }
        #   
        #   if (any(vektor == 0)) {
        #     break
        #   } #if zero in vector, ends - cause zero coverage
        #   
        #   count2 <- count2 + 1
        #   if (mod_vektor >= avg_cov_peak)
        #   {
        #     trigger2 <- trigger2 + 1
        #   }
        # }
        
        # writing the extended CNV
        vektor1 <- ( start - count1 * step):(start - 1)
        vektor <- peaksCELL[[i, 1]][, 1]
        vektor2 <- (stop + 1):(stop + count2 * step)
        
        matrixtemp <- cbind(c(vektor1, vektor, vektor2), coverageSignal[c(vektor1, vektor, vektor2)])
        peaksPolished[[i, 1]] <- matrixtemp
        CNVcoverage[i] <- avg_cov_peak
        CNVcoverageRelative[i] <- avg_cov_peak/(averageCoverage/100)
      } else{
        #NOT DELETIONS /\
        #SLOPE UPWARDS peak start /
        count1 <- 0
        trigger1 <- 0
        while (trigger1 < TRIGGER) {
          if (start <= 1 | (start - count1 * step - step) <= 1)
          {
            break
          }
          else
          {
            vektor <- coverageSignal[(start - count1 * step - step) : (start - count1 * step)]
            # mod_vektor <- mode(vektor) #matlab
            # mod_vektor <- sort(table(vektor),decreasing=TRUE)[1]
          }
          if (any(vektor == 0)) {
            break
          } #if zero in vector, ends - cause zero coverage
          
          slope <- (vektor[1] - vektor[length(vektor)]) / -length(vektor)
          
          count1 <- count1 + 1
          if (slope <= 0) {
            trigger1 <- trigger1 + 1
          }
        }
        
        # SLOPE DOWNWARDS peak stop \
        count2 <- 0
        trigger2 <- 0
        while (trigger2 < TRIGGER) {
          if (stop >= length(coverageSignal) |
              (stop + count2 * step + step) >= length(coverageSignal))
          {
            break
          }
          else{
            vektor <- coverageSignal[(stop + count2 * step) : (stop + count2 * step + step)]
            
            # mod_vektor <- mode(vektor); #matlab
            mod_vektor <- sort(table(vektor), decreasing = TRUE)[1]
          }
          
          if (any(vektor == 0)) {
            break
          } #if zero in vector, ends - cause zero coverage
          
          slope <- (vektor[1] - vektor[length(vektor)]) / -length(vektor)
          
          count2 <- count2 + 1
          if (slope >= 0)
          {
            trigger2 <- trigger2 + 1
          }
        }
        
        # writing the extended CNV
        vektor1 <- (start - count1 * step) : (start - 1)
        vektor <- peaksCELL[[i, 1]][, 1]
        vektor2 <- (stop + 1) : (stop + count2 * step)
        
        matrixtemp <-
          cbind(c(vektor1, vektor, vektor2), coverageSignal[c(vektor1, vektor, vektor2)])
        peaksPolished[[i, 1]] <- matrixtemp
        CNVcoverage[i] <- avg_cov_peak
        CNVcoverageRelative[i] <- avg_cov_peak/(averageCoverage/100)
      }
    }
  }
  
  ################################################################################
  ## SAVE CNV EVENTS into dataframe
  # coverage$isCNV <- FALSE
  DUP_DF <- data.frame(ID=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), TYPE=character())
  
  if(length(peaksPolished)>0){
    for(i in 1:length(peaksPolished)){
      # coverage$isCNV[peaksPolished[[i, 1]][,1]] <- TRUE
      start <-  peaksPolished[[i, 1]] [1, 1]
      nrows <- nrow(peaksPolished[[i, 1]])
      stop <- peaksPolished[[i, 1]] [nrows, 1]
      
      DUP_DF[i, "ID"] <- i
      DUP_DF[i, "START"] <- start
      DUP_DF[i, "END"] <- stop
      DUP_DF[i, "LENGTH"] <- (stop - start)
      DUP_DF[i, "COVERAGE"] <- CNVcoverage[i] #/ CNVcoverageRelative[i]
      DUP_DF[i, "TYPE"] <- "DUP"
      
    }
    ## check coordinates START<END
    DUP_DF <- subset(DUP_DF, START <= END)
  }
  
  
  return(DUP_DF)
}

