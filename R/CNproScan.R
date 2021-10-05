#' CNproScan - CNV detection for bacteria
#' This function detects and annotate CNVs in bacterial genomes. It uses GESD outliers detection 
#' and detection of discordant read-pairs for annotation. 
#' Last Update: 5/10/2021
#' 
#' @importFrom parallel splitIndices makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel
#' @importFrom seqinr getLength read.fasta
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @param infile Path to the input file
#' @return A vector of outliers indices from input
#' @export
#' 
CNproScanCNV <- function(coveragefile,bamFile,fastaFile,cores=2){
  ################################################################################
  
  # CONSTANT PARAMETERS
  peakDistanceThreshold<-20
  step<-11
  TRIGGER<-5
  
  ################################################################################
  ##  IMPORT COVERAGE FILE
  coverage <- read.csv(coveragefile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  header <- c('ID', 'POS', 'COVERAGE') # only 3 columns there
  colnames(coverage) <- header
  coverage <- coverage[,-1] #keeps only POS and COVERAGE
  
  averageCoverage<-mean(coverage$COVERAGE)
  
  ################################################################################
  ##  READ FASTA FILE
  referenceLength<-seqinr::getLength(seqinr::read.fasta(fastaFile,seqonly = TRUE))
  
  ################################################################################
  ## DELETING crossovers at the start and at the end and replace them by mean of whole coverage
  cutoff <- 500
  coverage$COVERAGE[1:cutoff] <- rep(mean(coverage$COVERAGE), times=length(coverage$COVERAGE[1:cutoff] ))
  end <- length(coverage$COVERAGE)
  coverage$COVERAGE[(end-cutoff):end] <- rep(mean(coverage$COVERAGE), times=length(coverage$COVERAGE[(end-cutoff):end] ))
  
  ################################################################################
  ## ZERO COVERAGE -> DELETIONS, 
  deletions<-which(coverage$COVERAGE==0)
  
  ################################################################################
  ## MERGING CLOSE DELETIONS TOGETHER
  idxDEL<-sort(deletions)
  differences<-diff(idxDEL)
  kk<-which(differences>1 & differences<peakDistanceThreshold)
  
  if(length(kk)){
    for (i in seq(1,length(kk))){
      CNVdistance <- idxDEL[kk[i]+1] - idxDEL[kk[i]]
      
      vec <- seq(idxDEL[kk[i]]+1, idxDEL[kk[i]]+CNVdistance-1) # coordinates to fill in the gap between 2 outliers
      #concatenate in between
      idxDEL <-c ( idxDEL[1:kk[i]] ,vec, idxDEL[(kk[i]+1):length(idxDEL)] )
      
      korekce <- CNVdistance-1 #correction - idx longer by one
      
      kk <- kk+korekce # correction so for loop is not biased by increased kk
    }
  }
  
  ################################################################################
  ## PREPROCESSING of DELETIONS into cell datatype
  idx <- idxDEL
  outliers_value <- coverage$COVERAGE[idx]
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
    ## SAVE into dataframe
    DEL_DF <- data.frame(ID=as.character(),START=integer(),END=integer(),LENGTH=integer(),COVERAGE=integer(),TYPE=character())
    for(i in 1:length(peaks)){
      start <-  peaks[[i, 1]] [1, 1]
      nrows <- nrow(peaks[[i, 1]])
      stop <- peaks[[i, 1]] [nrows, 1]
      DEL_DF[i,"ID"] <- i
      DEL_DF[i,"START"] <- start
      DEL_DF[i,"END"] <- stop
      DEL_DF[i,"LENGTH"] <- (stop - start)
      DEL_DF[i,"COVERAGE"] <- mean(coverage$COVERAGE[start:stop])
      DEL_DF[i,"TYPE"] <- "DEL"
    }
    ## DELETIONS COVERAGE VALUE REPLACED BY MEAN - INFLUENCES GESD OUTLIERS DETECTION
    coverage$COVERAGE[idxDEL] <- mean(coverage$COVERAGE)
  }
  # if only single deletion is detected
  if(length(k)==0 & length(idxDEL)!=0){
    DEL_DF <- data.frame(ID=as.character(),START=integer(),END=integer(),LENGTH=integer(),COVERAGE=integer(),TYPE=character())
    start <-  idxDEL[1]
    stop <- idxDEL[length(idxDEL)]
    DEL_DF[i,"ID"] <- 1
    DEL_DF[i,"START"] <- start
    DEL_DF[i,"END"] <- stop
    DEL_DF[i,"LENGTH"] <- (stop - start)
    DEL_DF[i,"COVERAGE"] <- mean(coverage$COVERAGE[start:stop])
    DEL_DF[i,"TYPE"] <- "DEL"
  }
  ################################################################################
  ## CNV EVENTS
  ## PARALLEL LOOP OUTLIERS DETECTION
  
  # cores <- parallel::detectCores(logical = FALSE) #should be odd?
  chunks <- parallel::splitIndices(length(coverage$COVERAGE), ncl = cores)
  
  #create the cluster
  my.cluster <- parallel::makeCluster(cores, type = "PSOCK")
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #  .packages = c("CNproScan"),
  Outliers <- foreach::foreach(i=1:cores, .export=c("gesd","mzscore"), .combine='c', .inorder=FALSE, .multicombine=FALSE) %dopar% {
    VECTOR<-coverage[chunks[[i]],]
    estimatedOutliers <- mzscore(VECTOR$COVERAGE)
    # if estimatedOutliers valid run GESD, check if is valid, else return empty vector
    if (!length(estimatedOutliers) == 0)
    {
      candidateOutliers <- gesd(VECTOR, length(estimatedOutliers))
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
  
  Outliers <- Outliers[Outliers!=0] #remove zeros
  
  ##############################################################################
  ## MERGING CLOSE CNV EVENTS TOGETHER
  idxGESD<-sort(Outliers)
  differences<-diff(idxGESD)
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
  outliers_value <- coverage$COVERAGE[idx]
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
  
  ################################################################################
  ## CNV EVENTS BORDERS DETECTION based on slope of a peak
  # works only on POSITIVE PEAKS, NOT NEGATIVE PEAKS - DELETIONS ?
  peaksPolished <- vector("list", length = length(peaks))
  dim(peaksPolished) <- c(length(peaks), 1)
  
  CNVcoverage <- c()
  CNVcoverageRelative <- c()
  
  coverageSignal <- coverage$COVERAGE
  
  for (i in seq(1, length(peaks))) {
    start <-  peaks[[i, 1]] [1, 1]
    nrows <- nrow(peaks[[i, 1]])
    stop <- peaks[[i, 1]] [nrows, 1]
    
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
      vektor1 <-( start - count1 * step):(start - 1)
      vektor <- peaks[[i, 1]][, 1]
      vektor2 <- (stop + 1):(stop + count2 * step)
      
      matrixtemp <-
        cbind(c(vektor1, vektor, vektor2), coverageSignal[c(vektor1, vektor, vektor2)])
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
      vektor <- peaks[[i, 1]][, 1]
      vektor2 <- (stop + 1) : (stop + count2 * step)
      
      matrixtemp <-
        cbind(c(vektor1, vektor, vektor2), coverageSignal[c(vektor1, vektor, vektor2)])
      peaksPolished[[i, 1]] <- matrixtemp
      CNVcoverage[i] <- avg_cov_peak
      CNVcoverageRelative[i] <- avg_cov_peak/(averageCoverage/100)
    }
  }
  
  ################################################################################
  ## SAVE CNV EVENTS into dataframe
  # coverage$isCNV <- FALSE
  DUP_DF <- data.frame(ID=as.character(),START=integer(),END=integer(),LENGTH=integer(),COVERAGE=integer(),TYPE=character())
  
  for(i in 1:length(peaksPolished)){
    # coverage$isCNV[peaksPolished[[i, 1]][,1]] <- TRUE
    start <-  peaksPolished[[i, 1]] [1, 1]
    nrows <- nrow(peaksPolished[[i, 1]])
    stop <- peaksPolished[[i, 1]] [nrows, 1]
    
    DUP_DF[i,"ID"] <- i
    DUP_DF[i,"START"] <- start
    DUP_DF[i,"END"] <- stop
    DUP_DF[i,"LENGTH"] <- (stop - start)
    DUP_DF[i,"COVERAGE"] <- CNVcoverage[i] #/ CNVcoverageRelative[i]
    DUP_DF[i,"TYPE"] <- "DUP"
    
  }
  
  ## check coordinates START<END
  DUP_DF <- subset(DUP_DF, START<=END)
  ## MERGE DELETIONS AND DUPLICATIONS
  CNV_DF <- rbind(DEL_DF,DUP_DF)
  # DELETE SINGLE BASE EVENTS
  CNV_DF <- subset(CNV_DF, LENGTH>=1)
  # CHECK FOR CNV < averageCoverage - those are deletions
  CNV_DF[CNV_DF$COVERAGE<averageCoverage,"TYPE"]<-"DEL"
  # R GARBAGE COLLECTOR
  gc()
  
  ################################################################################
  ## GET BASIC BAM INFO
  flag<-Rsamtools::scanBamFlag(isPaired = TRUE,isUnmappedQuery=FALSE,isProperPair=TRUE)
  parameters <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize","rname"),flag=flag)
  bam <- Rsamtools::scanBam(bamFile, param=parameters)[[1]]
  dfBAM<-as.data.frame(bam) #create data.frame
  
  ##GET SEQ HEADER aka REFERENCE 
  refHeader<-as.character(dfBAM$rname[1])
  ## DELETING crossovers at the start and at the end 
  cutoff <- 500 #defined above too
  dfBAM <- dfBAM[dfBAM$pos>cutoff,]
  dfBAM <- dfBAM[dfBAM$pos<(referenceLength-cutoff),]
  ## DELETE INSER SIZE OUTLIERS >10000
  isizecutoff <- 10000
  dfBAM <- dfBAM[abs(dfBAM$isize)<isizecutoff,]
  ## INSERT SIZE/TLEN - isize
  insertSize <- median(abs(dfBAM$isize))
  insertSize_1Q <- summary(abs(dfBAM$isize))[[2]]
  insertSize_3Q <- summary(abs(dfBAM$isize))[[5]]
  insertSize_IQR <- IQR(abs(dfBAM$isize))
  ## READ LENGTH - qwidth: The width of the query, as calculated from the cigar encoding; normally equal to the width of the query returned in seq.
  readLength <- median(abs(dfBAM$qwidth))
  # readLengthSD <- sd(abs(dfBAM$qwidth))
  
  ################################################################################
  ## SCAN BAM in CNV REGIONS
  
  for(i in 1:nrow(CNV_DF)){
    start <- CNV_DF[i,"START"]
    stop <- CNV_DF[i,"END"]
    region <- GenomicRanges::GRanges(refHeader, IRanges(start-insertSize, stop+insertSize))
    
    #-----------------------------------------------------------------------------
    ##  ACCESS BAM - TOTAL READS COUNT
    flagTOTAL <- Rsamtools::scanBamFlag(isPaired = TRUE)
    parametersTOTAL <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize"),which=region, flag=flagTOTAL)
    bamTOTAL <- Rsamtools::scanBam(bamFile, param=parametersTOTAL)[[1]]
    dfBAMTOTAL <- as.data.frame(bamTOTAL)
    dfBAMTOTAL <- na.omit(dfBAMTOTAL)
    if(!is.null(nrow(dfBAMTOTAL))){reads_TOTAL <- nrow(dfBAMTOTAL)}else{reads_TOTAL <- 0}
    
    #-----------------------------------------------------------------------------
    ##  ACCESS BAM - DELETIONS + INSERT SIZE OUTLIERS
    # flagDel<-scanBamFlag(isPaired = TRUE,isProperPair=TRUE)
    # flagDel<-scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE) # reads orientation -/+
    flagDel <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isFirstMateRead=TRUE,isMateMinusStrand=TRUE) # reads orientation +/-
    parametersDel <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize"),which=region, flag=flagDel)
    bamDel <- Rsamtools::scanBam(bamFile, param=parametersDel)[[1]]
    dfBAMDEL <- as.data.frame(bamDel)
    dfBAMDEL <- na.omit(dfBAMDEL)
    # if(!is.null(nrow(dfBAMDEL))){reads_DEL <- nrow(dfBAMDEL)}else{reads_DEL <- 0}
    
    
    # %% Circular genome correction
    if(!is.null(nrow(dfBAMDEL))){
      # % positive values
      REF <- referenceLength/2
      dfBAMDEL[dfBAMDEL$isize>REF,"isize"] <- dfBAMDEL[dfBAMDEL$isize>REF,"isize"] - 2*dfBAMDEL[dfBAMDEL$isize>REF,"qwidth"] - referenceLength
      # % negative values
      REF <- referenceLength/-2
      dfBAMDEL[dfBAMDEL$isize<REF,"isize"] <- dfBAMDEL[dfBAMDEL$isize<REF,"isize"] + 2*dfBAMDEL[dfBAMDEL$isize<REF,"qwidth"] + referenceLength
      #find outliers in TLEN fragment size 1.5xIQR rule
      out_ind <- which(abs(dfBAMDEL$isize)>(insertSize_3Q+1.5*insertSize_IQR))
      dfBAMDEL_INCREASE<-dfBAMDEL[out_ind,]
      out_ind <- which(abs(dfBAMDEL$isize)<(insertSize_1Q-1.5*insertSize_IQR))
      dfBAMDEL_DECREASE<-dfBAMDEL[out_ind,]
    }
    if(!is.null(nrow(dfBAMDEL_INCREASE))){reads_count_INCREASE <- nrow(dfBAMDEL_INCREASE)}else{reads_count_INCREASE <- 0}
    if(!is.null(nrow(dfBAMDEL_DECREASE))){reads_count_DECREASE <- nrow(dfBAMDEL_DECREASE)}else{reads_count_DECREASE <- 0}
    
    #-----------------------------------------------------------------------------
    ## ACCESS BAM TANDEM DUPLICATIONS - reads orientation -/+
    flagTD1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
    parametersTD1 <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","pos","qwidth","mpos"), which=region, flag=flagTD1)
    bamTD1 <- Rsamtools::scanBam(bamFile, param=parametersTD1)[[1]]
    dfbamTD1<-as.data.frame(bamTD1)
    dfbamTD1 <- na.omit(dfbamTD1)
    if(!is.null(nrow(dfbamTD1))){reads_count_TD <- nrow(dfbamTD1)}else{reads_count_TD <- 0}
    
    
    #-----------------------------------------------------------------------------
    # ACCESS BAM INTERSPERSED DUPLICATIONS INVERSED - reads orientation -/- OR +/+  AND INCREASED isize
    # reads orientation -/-
    flagINTERDUP_INV1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isMateMinusStrand=TRUE)
    parametersINTERDUP_INV1 <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","pos","qwidth","mpos","isize"), which=region, flag=flagINTERDUP_INV1)
    bamINTERDUP_INV1<- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_INV1)[[1]]
    dfbamINTERDUP_INV1<-as.data.frame(bamINTERDUP_INV1)
    dfbamINTERDUP_INV1 <- na.omit(dfbamINTERDUP_INV1)
    
    # %% Circular genome correction
    if(!is.null(nrow(dfbamINTERDUP_INV1))){
      # % positive values
      REF <- referenceLength/2
      dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize>REF,"isize"] <- dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize>REF,"isize"] - 2*dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize>REF,"qwidth"] - referenceLength
      # % negative values
      REF <- referenceLength/-2
      dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize<REF,"isize"] <- dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize<REF,"isize"] + 2*dfbamINTERDUP_INV1[dfbamINTERDUP_INV1$isize<REF,"qwidth"] + referenceLength
      #find outliers in TLEN fragment size #>insertSize_3Q
      out_ind <- which(abs(dfbamINTERDUP_INV1$isize)>(insertSize_3Q+1.5*insertSize_IQR))
      dfbamINTERDUP_INV1<-dfbamINTERDUP_INV1[out_ind,]
    }
    
    if(!is.null(nrow(dfbamINTERDUP_INV1))){reads_count_INTERDUP_INV1 <- nrow(dfbamINTERDUP_INV1)}else{reads_count_INTERDUP_INV1 <- 0}
    
    # reads orientation +/+
    flagINTERDUP_INV2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isMateMinusStrand=FALSE)
    parametersINTERDUP_INV2 <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","pos","qwidth","mpos","isize"),which=region,flag=flagINTERDUP_INV2)
    bamINTERDUP_INV2 <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_INV2)[[1]]
    dfbamINTERDUP_INV2 <- as.data.frame(bamINTERDUP_INV2)
    dfbamINTERDUP_INV2 <- na.omit(bamINTERDUP_INV2)
    
    # %% Circular genome correction
    if(!is.null(nrow(dfbamINTERDUP_INV2))){
      # % positive values
      REF <- referenceLength/2
      bamINTERDUP_INV2[bamINTERDUP_INV2$isize>REF,"isize"] <- bamINTERDUP_INV2[bamINTERDUP_INV2$isize>REF,"isize"] - 2*bamINTERDUP_INV2[bamINTERDUP_INV2$isize>REF,"qwidth"] - referenceLength
      # % negative values
      REF <- referenceLength/-2
      bamINTERDUP_INV2[bamINTERDUP_INV2$isize<REF,"isize"] <- bamINTERDUP_INV2[bamINTERDUP_INV2$isize<REF,"isize"] + 2*bamINTERDUP_INV2[bamINTERDUP_INV2$isize<REF,"qwidth"] + referenceLength
      #find outliers in TLEN fragment size #>insertSize_3Q
      out_ind <- which(abs(bamINTERDUP_INV2$isize)>(insertSize_3Q+1.5*insertSize_IQR))
      bamINTERDUP_INV2<-bamINTERDUP_INV2[out_ind,]
    }
    
    if(!is.null(nrow(dfbamINTERDUP_INV2))){reads_count_INTERDUP_INV2 <- nrow(dfbamINTERDUP_INV2)}else{reads_count_INTERDUP_INV2 <- 0}
    
    reads_count_INTERDUP_INV <- reads_count_INTERDUP_INV2 + reads_count_INTERDUP_INV1
    
    #-----------------------------------------------------------------------------
    # ACCESS BAM INTERSPERSED DUPLICATIONS DIRECT - reads orientation +/- OR -/+  AND INCREASED isize
    # reads orientation +/-
    flagINTERDUP_DIR1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isFirstMateRead=TRUE,isMateMinusStrand=TRUE)
    parametersINTERDUP_DIR1 <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","pos","qwidth","mpos","isize"), which=region, flag=flagINTERDUP_DIR1)
    bamINTERDUP_DIR1 <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_DIR1)[[1]]
    dfbamINTERDUP_DIR1 <- as.data.frame(bamINTERDUP_DIR1)
    dfbamINTERDUP_DIR1 <- na.omit(dfbamINTERDUP_DIR1)
    
    # %% Circular genome correction
    if(!is.null(nrow(dfbamINTERDUP_DIR1))){
      # % positive values
      REF <- referenceLength/2
      dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize>REF,"isize"] <- dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize>REF,"isize"] - 2*dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize>REF,"qwidth"] - referenceLength
      # % negative values
      REF <- referenceLength/-2
      dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize<REF,"isize"] <- dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize<REF,"isize"] + 2*dfbamINTERDUP_DIR1[dfbamINTERDUP_DIR1$isize<REF,"qwidth"] + referenceLength
      #find outliers in TLEN fragment size #>insertSize_3Q
      out_ind <- which(abs(dfbamINTERDUP_DIR1$isize)>(insertSize_3Q+1.5*insertSize_IQR))
      dfbamINTERDUP_DIR1<-dfbamINTERDUP_DIR1[out_ind,]
    }
    
    if(!is.null(nrow(dfbamINTERDUP_DIR1))){reads_count_INTERDUP_DIR1 <- nrow(dfbamINTERDUP_DIR1)}else{reads_count_INTERDUP_DIR1 <- 0}
    
    # reads orientation -/+
    flagINTERDUP_DIR2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
    parametersINTERDUP_DIR2 <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","pos","qwidth","mpos","isize"),which=region,flag=flagINTERDUP_DIR2)
    bamINTERDUP_DIR2<- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_DIR2)[[1]]
    dfbamINTERDUP_DIR2<-as.data.frame(bamINTERDUP_DIR2)
    dfbamINTERDUP_DIR2 <- na.omit(dfbamINTERDUP_DIR2)
    
    # %% Circular genome correction
    if(!is.null(nrow(dfbamINTERDUP_DIR2))){
      # % positive values
      REF <- referenceLength/2
      dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize>REF,"isize"] <- dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize>REF,"isize"] - 2*dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize>REF,"qwidth"] - referenceLength
      # % negative values
      REF <- referenceLength/-2
      dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize<REF,"isize"] <- dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize<REF,"isize"] + 2*dfbamINTERDUP_DIR2[dfbamINTERDUP_DIR2$isize<REF,"qwidth"] + referenceLength
      #find outliers in TLEN fragment size #>insertSize_3Q
      out_ind <- which(abs(dfbamINTERDUP_DIR2$isize)>(insertSize_3Q+1.5*insertSize_IQR))
      dfbamINTERDUP_DIR2<-dfbamINTERDUP_DIR2[out_ind,]
    }
    
    if(!is.null(nrow(dfbamINTERDUP_DIR2))){reads_count_INTERDUP_DIR2 <- nrow(dfbamINTERDUP_DIR2)}else{reads_count_INTERDUP_DIR2 <- 0}
    
    reads_count_INTERDUP_DIR <- reads_count_INTERDUP_DIR2 + reads_count_INTERDUP_DIR1
    
    # #print absolute numbers
    # CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
    # CNV_DF[i,"reads_INSERT_INCREASE"] <- reads_count_INCREASE
    # CNV_DF[i,"reads_INSERT_DECREASE"] <- reads_count_DECREASE
    # CNV_DF[i,"reads_count_TAD_DUP"] <- reads_count_TD
    # CNV_DF[i,"reads_count_INTER_DUP_DIR"] <- reads_count_INTERDUP_DIR
    # CNV_DF[i,"reads_count_INTER_DUP_INV"] <- reads_count_INTERDUP_INV
    
    #print percents
    CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
    CNV_DF[i,"reads_INSERT_INCREASE"] <- reads_count_INCREASE/(reads_TOTAL/100)
    CNV_DF[i,"reads_INSERT_DECREASE"] <- reads_count_DECREASE/(reads_TOTAL/100)
    CNV_DF[i,"reads_count_TAD_DUP"] <- reads_count_TD/(reads_TOTAL/100)
    CNV_DF[i,"reads_count_INTER_DUP_DIR"] <- reads_count_INTERDUP_DIR/(reads_TOTAL/100)
    CNV_DF[i,"reads_count_INTER_DUP_INV"] <- reads_count_INTERDUP_INV/(reads_TOTAL/100)
    
    #TYPE SELECTION
    if(CNV_DF[i,"TYPE"] == "DEL"){CNV_DF[i,"SUBTYPE"] <- "DEL"}
    if(CNV_DF[i,"TYPE"] == "DUP") {CNV_DF[i,"SUBTYPE"] <- "TANDEM"}
    if(CNV_DF[i,"TYPE"] == "DUP" &  CNV_DF[i,"reads_count_INTER_DUP_DIR"]>CNV_DF[i,"reads_count_TAD_DUP"] /4) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_DIRECT"}
    if(CNV_DF[i,"TYPE"] == "DUP" &  CNV_DF[i,"reads_count_INTER_DUP_INV"]>CNV_DF[i,"reads_count_TAD_DUP"] /4) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_INVERSED"}
    
  }
  
  return(CNV_DF)
  
}
