#' cnvDiscordantReads
#' Goes over BAM file and search for discordant reads. Returns input dataframe with added information. 
#' 
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @param bamFile path to the BAM file, sorted, indexed
#' @param CNV_DF dataframe of already detected CNVs
#' @param refHeader name of chromosome or conrig
#' @param referenceLength lengh of given contig/chromosome
#' @return CNV dataframe
#' @export
#' @noRd
#' 
cnvDiscordantReads <- function(bamFile, CNV_DF, refHeader, referenceLength){
  
  # DISCORDANT READS
  ## GET BASIC BAM INFO
  flag<-Rsamtools::scanBamFlag(isPaired = TRUE,isUnmappedQuery=FALSE,isProperPair=TRUE)
  parameters <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize","rname"),flag=flag)
  bam <- Rsamtools::scanBam(bamFile, param=parameters)[[1]]
  dfBAM<-as.data.frame(bam) #create data.frame
  
  ##GET ONLY SELECTED CONTIG 
  dfBAM <- dfBAM[dfBAM$rname== refHeader,]
  
  ## DELETING crossovers at the start and at the end 
  cutoff <- 100 #defined above too
  dfBAM <- dfBAM[dfBAM$pos>cutoff,]
  dfBAM <- dfBAM[dfBAM$pos<(referenceLength-cutoff),]
  ## DELETE INSER SIZE OUTLIERS >10000
  # isizecutoff <- 10000
  # dfBAM <- dfBAM[abs(dfBAM$isize)<isizecutoff,]
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
  
  if(nrow(CNV_DF)>0){
    for(i in 1:nrow(CNV_DF)){
      start <- CNV_DF[i,"START"]
      stop <- CNV_DF[i,"END"]
      region <- GenomicRanges::GRanges(refHeader, IRanges::IRanges(start-insertSize, stop+insertSize))
      
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
      CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
      CNV_DF[i,"reads_INSERT_INCREASE"] <- reads_count_INCREASE
      CNV_DF[i,"reads_INSERT_DECREASE"] <- reads_count_DECREASE
      CNV_DF[i,"reads_count_TAD_DUP"] <- reads_count_TD
      CNV_DF[i,"reads_count_INTER_DUP_DIR"] <- reads_count_INTERDUP_DIR
      CNV_DF[i,"reads_count_INTER_DUP_INV"] <- reads_count_INTERDUP_INV
      
      #print Percentages
      # CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
      # CNV_DF[i,"reads_INSERT_INCREASE"] <- reads_count_INCREASE/(reads_TOTAL/100)
      # CNV_DF[i,"reads_INSERT_DECREASE"] <- reads_count_DECREASE/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_TAD_DUP"] <- reads_count_TD/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_INTER_DUP_DIR"] <- reads_count_INTERDUP_DIR/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_INTER_DUP_INV"] <- reads_count_INTERDUP_INV/(reads_TOTAL/100)
      
      #TYPE SELECTION
      if(CNV_DF[i,"TYPE"] == "DEL"){CNV_DF[i,"SUBTYPE"] <- "DEL"}
      if(CNV_DF[i,"TYPE"] == "DUP") {CNV_DF[i,"SUBTYPE"] <- "TANDEM"}
      if(CNV_DF[i,"TYPE"] == "DUP" &  CNV_DF[i,"reads_count_INTER_DUP_DIR"]>CNV_DF[i,"reads_count_TAD_DUP"] /4) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_DIRECT"}
      if(CNV_DF[i,"TYPE"] == "DUP" &  CNV_DF[i,"reads_count_INTER_DUP_INV"]>CNV_DF[i,"reads_count_TAD_DUP"] /4) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_INVERSED"}
      
    }
  } else{CNV_DF <- data.frame(ID=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), TYPE=character(),
                              reads_TOTAL=integer(),reads_INSERT_INCREASE=integer(),reads_INSERT_DECREASE=integer(),reads_count_TAD_DUP=integer(),
                              reads_count_INTER_DUP_DIR=integer(),reads_count_INTER_DUP_INV=integer(),SUBTYPE=as.character()
  )
  
        CNV_DF[1,] <- NA
  }
  
  return(CNV_DF)
}

