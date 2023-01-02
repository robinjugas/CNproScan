#' cnvDiscordantReads v2
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
  
  ## INSERT SIZE/TLEN - isize
  insertSize <- median(abs(dfBAM$isize))
  insertSize_1Q <- summary(abs(dfBAM$isize))[[2]]
  insertSize_3Q <- summary(abs(dfBAM$isize))[[5]]
  insertSize_IQR <- IQR(abs(dfBAM$isize))
  ## READ LENGTH - qwidth: The width of the query, as calculated from the cigar encoding; normally equal to the width of the query returned in seq.
  readLength <- median(abs(dfBAM$qwidth))
  
  ################################################################################
  ## SCAN BAM in CNV REGIONS
  
  if(nrow(CNV_DF)>0){
    for(i in 1:nrow(CNV_DF)){
      
      start <- CNV_DF[i,"START"]
      stop <- CNV_DF[i,"END"]
      # Granges boundaries 
      STARTGRANGES <- start-insertSize
      STOPGRANGES <- stop+insertSize
      if (STARTGRANGES<1){STARTGRANGES <- 1}
      if(STOPGRANGES>referenceLength){STOPGRANGES <- referenceLength}
      region <- GenomicRanges::GRanges(refHeader, IRanges::IRanges(start-insertSize, stop+insertSize))
      
      CNVLENGTH <- CNV_DF[i,"LENGTH"]
      #-----------------------------------------------------------------------------
      ##  ACCESS BAM - TOTAL READS COUNT
      flagTOTAL <- Rsamtools::scanBamFlag(isPaired = TRUE)
      parametersTOTAL <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize"),which=region, flag=flagTOTAL)
      bamTOTAL <- Rsamtools::scanBam(bamFile, param=parametersTOTAL)[[1]]
      dfBAMTOTAL <- as.data.frame(bamTOTAL)
      dfBAMTOTAL <- na.omit(dfBAMTOTAL)
      if(!is.null(nrow(dfBAMTOTAL))){reads_TOTAL <- nrow(dfBAMTOTAL)}else{reads_TOTAL <- 0}
      
      ################################################################################################################################################################
      #-----------------------------------------------------------------------------
      ##  ACCESS BAM - DELETIONS - reads orientation +/- -/+  AND higher INSERT SIZE
      # reads orientation  +/-
      flagDel1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isFirstMateRead=TRUE,isMateMinusStrand=TRUE) # reads orientation +/-
      parametersDel1 <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize"),which=region, flag=flagDel1)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersDel1)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+1.5*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_INCREASE))){reads_supporting_DELETION1 <- nrow(dfBAM_INCREASE)}else{reads_supporting_DELETION1 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_INCREASE))){reads_supporting_DELETION1 <- reads_supporting_DELETION1+nrow(dfBAM_ASCNVLENGTH)}
      
      # reads orientation  -/+
      flagDel2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
      parametersDel2 <- Rsamtools::ScanBamParam(what=c("pos", "qwidth","isize"),which=region, flag=flagDel2)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersDel2)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+1.5*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_INCREASE))){reads_supporting_DELETION2 <- nrow(dfBAM_INCREASE)}else{reads_supporting_DELETION2 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_INCREASE))){reads_supporting_DELETION2 <- reads_supporting_DELETION2+nrow(dfBAM_ASCNVLENGTH)}
      
      reads_supporting_DELETION <- reads_supporting_DELETION1 + reads_supporting_DELETION2
      
      ################################################################################################################################################################
      #-----------------------------------------------------------------------------
      ## ACCESS BAM TANDEM DUPLICATIONS DIRECT - reads orientation -/+ or +/-  AND lower INSERT SIZE OR INCREASE AROUND SIZE of CNV length +- insertSize 
      # reads orientation  -/+
      flagTD_DIR1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
      parametersTD_DIR1 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagTD_DIR1)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersTD_DIR1)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)<(insertSize_1Q-1.5*insertSize_IQR))
        dfBAM_DECREASE<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADDIR1 <- nrow(dfBAM_DECREASE)}else{reads_supporting_TADDIR1 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADDIR1 <- reads_supporting_TADDIR1+nrow(dfBAM_ASCNVLENGTH)}
      
      
      # reads orientation  +/-
      flagTD_DIR2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
      parametersTD_DIR2 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagTD_DIR2)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersTD_DIR2)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)<(insertSize_1Q-1.5*insertSize_IQR))
        dfBAM_DECREASE<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADDIR2 <- nrow(dfBAM_DECREASE)}else{reads_supporting_TADDIR2 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADDIR2 <- reads_supporting_TADDIR2+nrow(dfBAM_ASCNVLENGTH)}
      reads_supporting_TANDEM_DIRECT <- reads_supporting_TADDIR1 + reads_supporting_TADDIR2
      
      
      #-----------------------------------------------------------------------------
      ## ACCESS BAM TANDEM DUPLICATIONS INDIRECT - reads orientation -/- or +/+  AND lower INSERT SIZE
      # reads orientation  -/-
      flagTD_INDIR1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isMateMinusStrand=TRUE)
      parametersTD_INDIR1 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagTD_INDIR1)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersTD_INDIR1)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)<(insertSize_1Q-1.5*insertSize_IQR))
        dfBAM_DECREASE<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADINDIR1 <- nrow(dfBAM_DECREASE)}else{reads_supporting_TADINDIR1 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADINDIR1 <- reads_supporting_TADINDIR1+nrow(dfBAM_ASCNVLENGTH)}
      
      
      # reads orientation  +/+
      flagTD_INDIR2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isMateMinusStrand=FALSE)
      parametersTD_INDIR2 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagTD_INDIR2)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersTD_INDIR2)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule or as same as CNV length
        out_ind <- which(abs(dfBAM$isize)>(CNVLENGTH-2*insertSize) &  abs(dfBAM$isize)<(CNVLENGTH+2*insertSize))
        dfBAM_ASCNVLENGTH<-dfBAM[out_ind,]
        out_ind <- which(abs(dfBAM$isize)<(insertSize_1Q-1.5*insertSize_IQR))
        dfBAM_DECREASE<-dfBAM[out_ind,]
      }
      if(!is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADINDIR2 <- nrow(dfBAM_DECREASE)}else{reads_supporting_TADINDIR2 <- 0}
      if(!is.null(nrow(dfBAM_ASCNVLENGTH)) & !is.null(nrow(dfBAM_DECREASE))){reads_supporting_TADINDIR2 <- reads_supporting_TADINDIR2+nrow(dfBAM_ASCNVLENGTH)}
      
      reads_supporting_TANDEM_INDIRECT <- reads_supporting_TADINDIR1 + reads_supporting_TADINDIR2
      
      ################################################################################################################################################################
      #-----------------------------------------------------------------------------
      # ACCESS BAM INTERSPERSED DUPLICATIONS INVERSED - reads orientation -/- OR +/+  AND much higher INSERT SIZE
      # reads orientation -/-
      flagINTERDUP_INV1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isMateMinusStrand=TRUE)
      parametersINTERDUP_INV1 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagINTERDUP_INV1)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_INV1)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+10*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
      }
      
      if(!is.null(nrow(dfBAM_INCREASE))){reads_count_INTERDUP_INV1 <- nrow(dfBAM_INCREASE)}else{reads_count_INTERDUP_INV1 <- 0}
      
      # reads orientation +/+
      flagINTERDUP_INV2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isMateMinusStrand=FALSE)
      parametersINTERDUP_INV2 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"),which=region,flag=flagINTERDUP_INV2)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_INV2)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+10*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
      }
      
      if(!is.null(nrow(dfBAM_INCREASE))){reads_count_INTERDUP_INV2 <- nrow(dfBAM_INCREASE)}else{reads_count_INTERDUP_INV2 <- 0}
      
      reads_supporting_INTERSPERSED_INDIRECT <- reads_count_INTERDUP_INV2 + reads_count_INTERDUP_INV1
      
      #-----------------------------------------------------------------------------
      # ACCESS BAM INTERSPERSED DUPLICATIONS DIRECT - reads orientation +/- OR -/+  AND much higher INSERT SIZE
      # reads orientation +/-
      flagINTERDUP_DIR1 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=FALSE, isFirstMateRead=TRUE,isMateMinusStrand=TRUE)
      parametersINTERDUP_DIR1 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"), which=region, flag=flagINTERDUP_DIR1)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_DIR1)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+10*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
      }
      
      if(!is.null(nrow(dfBAM_INCREASE))){reads_count_INTERDUP_DIR1 <- nrow(dfBAM_INCREASE)}else{reads_count_INTERDUP_DIR1 <- 0}
      
      # reads orientation -/+
      flagINTERDUP_DIR2 <- Rsamtools::scanBamFlag(isPaired = TRUE, isMinusStrand=TRUE, isFirstMateRead=TRUE,isMateMinusStrand=FALSE)
      parametersINTERDUP_DIR2 <- Rsamtools::ScanBamParam(what = c("pos", "qwidth","isize"),which=region,flag=flagINTERDUP_DIR2)
      dfBAM <- Rsamtools::scanBam(bamFile, param=parametersINTERDUP_DIR2)[[1]]
      dfBAM <- as.data.frame(dfBAM)
      dfBAM <- na.omit(dfBAM)
      # %% Circular genome correction
      if(!is.null(nrow(dfBAM))){
        # % positive values
        REF <- referenceLength/2
        dfBAM[dfBAM$isize>REF,"isize"] <- dfBAM[dfBAM$isize>REF,"isize"] - 2*dfBAM[dfBAM$isize>REF,"qwidth"] - referenceLength
        # % negative values
        REF <- referenceLength/-2
        dfBAM[dfBAM$isize<REF,"isize"] <- dfBAM[dfBAM$isize<REF,"isize"] + 2*dfBAM[dfBAM$isize<REF,"qwidth"] + referenceLength
        #find outliers in TLEN fragment size 1.5xIQR rule
        out_ind <- which(abs(dfBAM$isize)>(insertSize_3Q+10*insertSize_IQR))
        dfBAM_INCREASE<-dfBAM[out_ind,]
      }
      
      if(!is.null(nrow(dfBAM_INCREASE))){reads_count_INTERDUP_DIR2 <- nrow(dfBAM_INCREASE)}else{reads_count_INTERDUP_DIR2 <- 0}
      
      reads_supporting_INTERSPERSED_DIRECT <- reads_count_INTERDUP_DIR2 + reads_count_INTERDUP_DIR1
      
      
      #-----------------------------------------------------------------------------
      # #print absolute numbers
      CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
      CNV_DF[i,"reads_supporting_Deletion"] <- reads_supporting_DELETION
      CNV_DF[i,"reads_supporting_Tandem_Direct"] <- reads_supporting_TANDEM_DIRECT
      CNV_DF[i,"reads_supporting_Tandem_Indirect"] <- reads_supporting_TANDEM_INDIRECT
      CNV_DF[i,"reads_supporting_Interspersed_Direct"] <- reads_supporting_INTERSPERSED_DIRECT
      CNV_DF[i,"reads_supporting_Interspersed_Indirect"] <- reads_supporting_INTERSPERSED_INDIRECT
      
      #print Percentages
      # CNV_DF[i,"reads_TOTAL"] <- reads_TOTAL
      # CNV_DF[i,"reads_INSERT_INCREASE"] <- reads_count_INCREASE/(reads_TOTAL/100)
      # CNV_DF[i,"reads_INSERT_DECREASE"] <- reads_count_DECREASE/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_TAD_DUP"] <- reads_count_TD/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_INTER_DUP_DIR"] <- reads_count_INTERDUP_DIR/(reads_TOTAL/100)
      # CNV_DF[i,"reads_count_INTER_DUP_INV"] <- reads_count_INTERDUP_INV/(reads_TOTAL/100)
      
      
      # SUBTYPE SELECTION
      signaturesSUM <- reads_supporting_DELETION+reads_supporting_TANDEM_DIRECT+reads_supporting_TANDEM_INDIRECT+reads_supporting_INTERSPERSED_DIRECT+reads_supporting_INTERSPERSED_INDIRECT
      tempX <- which.max(c(reads_supporting_TANDEM_DIRECT,reads_supporting_TANDEM_INDIRECT,reads_supporting_INTERSPERSED_DIRECT,reads_supporting_INTERSPERSED_INDIRECT))
      
      if(CNV_DF[i,"TYPE"] == "DEL") {CNV_DF[i,"SUBTYPE"] <- "DEL"}
      
      if(CNV_DF[i,"TYPE"] == "DUP" & tempX==1) {CNV_DF[i,"SUBTYPE"] <- "TANDEM_DIRECT"}
      if(CNV_DF[i,"TYPE"] == "DUP" & tempX==2) {CNV_DF[i,"SUBTYPE"] <- "TANDEM_INDIRECT"}
      if(CNV_DF[i,"TYPE"] == "DUP" & tempX==3) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_DIRECT"}
      if(CNV_DF[i,"TYPE"] == "DUP" & tempX==4) {CNV_DF[i,"SUBTYPE"] <- "INTERSPERSED_INDIRECT"}
      
      if(signaturesSUM < (reads_TOTAL/10) ) {CNV_DF[i,"SUBTYPE"] <- "FALSE POSITIVE"}
      
      #reorder
      CNV_DF <- CNV_DF[,c("ID",
                          "START",
                          "END",
                          "LENGTH",
                          "COVERAGE",
                          "TYPE",
                          "SUBTYPE",
                          "reads_TOTAL",
                          "reads_supporting_Deletion",
                          "reads_supporting_Tandem_Direct",
                          "reads_supporting_Tandem_Indirect",
                          "reads_supporting_Interspersed_Direct",
                          "reads_supporting_Interspersed_Indirect"
      )]
      
      
    }
  } else{CNV_DF <- data.frame(ID=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), TYPE=character(),SUBTYPE=as.character(),
                              reads_TOTAL=integer(),reads_supporting_Deletion=integer(),reads_supporting_Tandem_Direct=integer(),reads_supporting_Tandem_Indirect=integer(),
                              reads_supporting_Interspersed_Direct=integer(),reads_supporting_Interspersed_Indirect=integer())
  }
  
  
  
  return(CNV_DF)
}

