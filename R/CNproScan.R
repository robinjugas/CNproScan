#' CNproScan - CNV detection for bacteria genomes
#' This function detects and annotate CNVs in bacterial genomes. It uses GESD outliers detection 
#' and detection of discordant read-pairs for annotation. 
#' 
#' @importFrom seqinr getLength read.fasta
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @import data.table
#' @param coverageFile Path to the coverage file given from samtools depth. Zeroe values must be included, use -a switch "samtools depth -a".
#' @param bamFile Path to the BAM file, sorted and indexed. 
#' @param fastaFile Path to the reference FASTA file used to align sequencing reads. The headers should be the same through all input files (FASTA, BAM and COVERAGE files). 
#' @param GCnorm Should the GC normalization be performed. Logical value.  Default is TRUE. 
#' @param MAPnorm  Should the genome mappability normalization be performed. Default is FALSE. If TRUE, the argument bedgraphFile must be specified. 
#' @param bedgraphFile Path to the bedgraph file outputed from genmap tool. 
#' @param cores number of cores to be used by parallel package. Default is two. Decreases computation time with limited impact on RAM usage. 
#' @return A dataframe of detected CNVs and VCF file name as sample_cnproscan.vcf in the working directory. 
#' @export
#' 

CNproScanCNV <- function(coverageFile,bamFile,fastaFile,GCnorm=TRUE,MAPnorm=FALSE,ORICnorm=FALSE,bedgraphFile=NULL,oriCposition=0,cores=2){
  
  ################################################################################
  ## CHECK INPUTS
  stopifnot("`coverageFile` must be a character."=is.character(coverageFile))
  stopifnot("`bamFile` must be a character."=is.character(bamFile))
  stopifnot("`fastaFile` must be a character."=is.character(fastaFile))
  
  if (MAPnorm == TRUE & is.null(bedgraphFile)){
    stop("`bedgraphFile` must be specified if MAPnorm=TRUE.")
  }
  if(is.null(bedgraphFile)){bedgraphFile<-""}
  
  if (MAPnorm == TRUE & !file.exists(bedgraphFile)){
    stop("`bedgraphFile` file does not exist or is not a valid file.")
  }
  
  if (ORICnorm == TRUE & is.null(oriCposition)){
    stop("`oriCposition` has to be defined of oriC normalization is required.")
  }
  
  
  stopifnot("`coverageFile` file does not exist or is not a valid file."=file.exists(coverageFile))
  stopifnot("`bamFile` file does not exist or is not a valid file."=file.exists(bamFile))
  stopifnot("`fastaFile` file does not exist or is not a valid file."=file.exists(fastaFile))
  
  ################################################################################
  ## CONSTANT PARAMETERS
  peakDistanceThreshold <- 20 #distance between close peaks to be merged
  step <- 11 # number of bases to calculate slope
  TRIGGER <- 5 # number of changes of slope derivation
  
  ################################################################################
  ##  READ COVERAGE FILE
  coverage <- read.csv(coverageFile, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  if (ncol(coverage) != 3){
    stop("`coverageFile` must have 3 columns as outputed by 'samtools depth'. ")
  }
  header <- c('CHROM', 'POS', 'COVERAGE') # only 3 columns there
  colnames(coverage) <- header
  averageCoverage <- mean(coverage$COVERAGE)
  numberOfCoverageContigs <- length(unique(coverage$CHROM))
  CoverageContigList <- unique(coverage$CHROM)
  
  
  ################################################################################
  ##  READ FASTA FILE
  referenceLengthList <- seqinr::getLength(seqinr::read.fasta(fastaFile, seqonly=TRUE))
  numberOfFastaContigs <- length(referenceLengthList)
  fastaContiglist <- seqinr::getName(seqinr::read.fasta(fastaFile, seqonly=FALSE))
  
  if (length(referenceLengthList) == 0){
    stop("`fastaFile` must not be empty.")
  }
  
  ################################################################################
  ## READ BAM FILE
  flag <- Rsamtools::scanBamFlag(isPaired=TRUE, isUnmappedQuery=FALSE, isProperPair=TRUE)
  parameters <- Rsamtools::ScanBamParam(what=c("pos", "qwidth", "isize", "rname"), flag=flag)
  bam <- Rsamtools::scanBam(bamFile, param=parameters)[[1]]
  dfBAM <- as.data.frame(bam) #create data.frame
  numberOfBAMContigs <- length(unique(dfBAM$rname))
  BAMContigList <- as.character(unique(dfBAM$rname)) #GET CONTIGS
  rm(dfBAM, bam, parameters, flag)
  
  ################################################################################
  ## DIFFERENT CONTIGS NUMBER, also chromosome from coverage File may be numeric, while from BAM and FASTA an character
  
  # if not as string, make string
  # if(!is.character(CoverageContigList)){
  #   
  #   CoverageContigList <- as.character(CoverageContigList)
  #   }
  
  # CoverageContigList <- c(1,2,3,4,5)
  # CoverageContigList <- as.vector(CoverageContigList,"character")
  # CoverageContigList
  # 
  # CoverageContigList <- c(1,2,3,4,5)
  # CoverageContigList <- as.character(CoverageContigList)
  # CoverageContigList
  
  
  # keep only those in coverage File (aka CoverageContigList)
  l <- list(fastaContiglist,CoverageContigList,BAMContigList)
  # all(sapply(l, length) == length(l[[1]]))
  if(!all(sapply(l, length) == length(l[[1]]))){
    fastaContiglist <- intersect(fastaContiglist,CoverageContigList)
    BAMContigList <- intersect(BAMContigList,CoverageContigList)
  }
  
  CoverageContigList
  fastaContiglist
  BAMContigList
  
  numberOfFastaContigs <- length(fastaContiglist)
  numberOfBAMContigs <- length(BAMContigList)
  numberOfCoverageContigs <- length(CoverageContigList)
  
  
  ################################################################################
  ################################################################################
  ## RUN ONLY SINGLE CONTIG/CHROM
  if(numberOfFastaContigs == numberOfCoverageContigs & numberOfCoverageContigs == numberOfBAMContigs & numberOfCoverageContigs == 1){
    
    referenceLength <- as.numeric(referenceLengthList[1])
    refHeader <- BAMContigList[1] # get int
    fastaContig <- fastaContiglist[1]
    CoverageContig <- CoverageContigList[1]
    ## check for condition
    if(refHeader != fastaContig | fastaContig != CoverageContig | refHeader != CoverageContig){
      stop("There are different chromosome/contig names between BAM, FASTA and coverage file. Reference FASTA file used for alignment might be different than defined.")
    }
    
    ################################################################################
    ##  GC normalization
    if(GCnorm == TRUE){coverage$COVERAGE <- gcNormalization(coverage$COVERAGE, fastaFile, refHeader)}
    ################################################################################
    ##  MAPPABILITY normalization
    if(MAPnorm == TRUE & file.exists(bedgraphFile)){coverage$COVERAGE <- mappabilityNormalization(coverage$COVERAGE, bedgraphFile, refHeader)}
    ################################################################################
    ##  oriC normalization
    if(ORICnorm == TRUE & !is.null(oriCposition)){coverage$COVERAGE <- oricNormalization_v2(coverage$COVERAGE, referenceLength, oriCposition)}
    ################################################################################
    ## CNV EVENTS
    peaks <- cnvOutliers(coverage, cores, peakDistanceThreshold)
    ################################################################################
    ## CNV EVENTS BORDERS DETECTION based on slope of a peak
    CNV_DF <- cnvBoundaries(peaks, coverage, step, TRIGGER)
    ################################################################################
    ## MERGE
    CNV_DF <- subset(CNV_DF, LENGTH >= 1) # DELETE SINGLE BASE EVENTS 
    ################################################################################
    # DISCORDANT READS
    CNV_DF <- cnvDiscordantReads(bamFile, CNV_DF, refHeader, referenceLength)
    ################################################################################
    
    ## ADD CNV copy number
    if(nrow(CNV_DF)>0){
      ## ADD CONTIG NAME
      CNV_DF$CHROM <- refHeader
      CNV_DF$COPY_NUMBER <- round(CNV_DF$COVERAGE/averageCoverage)
    } else{CNV_DF <- data.frame(ID=as.character(), CHROM=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), COPY_NUMBER=integer(),TYPE=character(),SUBTYPE=as.character(),
                                reads_TOTAL=integer(),reads_supporting_Deletion=integer(),reads_supporting_Tandem_Direct=integer(),reads_supporting_Tandem_Indirect=integer(),
                                reads_supporting_Interspersed_Direct=integer(),reads_supporting_Interspersed_Indirect=integer())
    }
    
    #reorder columns
    if(nrow(CNV_DF)>0){
      CNV_DF <- CNV_DF[,c("ID",
                          "CHROM",
                          "START",
                          "END",
                          "LENGTH",
                          "COVERAGE",
                          "COPY_NUMBER",
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
    
    ## reformat and remove NAs CNVs
    CNV_DF <- CNV_DF[!(is.na(CNV_DF$START)),] # remove NA CNVs
    CNV_DF <- as.data.table(CNV_DF)
    CNV_DF <- unique(CNV_DF,by=c("CHROM","START","END","TYPE","SUBTYPE"))
    CNV_DF <- CNV_DF[order(CNV_DF$CHROM, CNV_DF$START),]
    
    
  }
  
  
  ################################################################################
  ################################################################################
  ## RUN MULTI - CONTIG/CHROM
  if(numberOfFastaContigs == numberOfCoverageContigs & numberOfCoverageContigs == numberOfBAMContigs & numberOfCoverageContigs > 1){   
    
    ## merge corresponding values in case they are differently ordered
    CONTIGS <- as.data.frame(cbind(fastaContiglist, referenceLengthList))
    CONTIGS <- merge(CONTIGS, as.data.frame(CoverageContigList), by.x="fastaContiglist", by.y="CoverageContigList") 
    CONTIGS <- merge(CONTIGS, as.data.frame(BAMContigList), by.x="fastaContiglist", by.y="BAMContigList") 
    colnames(CONTIGS) <- c("refHeader", "referenceLength")
    MERGED_CNV_DF <- data.frame()     # empty DF to initiate
    
    ## go over contigs/chromosome
    for (i in 1:nrow(CONTIGS)){
      
      referenceLength <- as.numeric(CONTIGS$referenceLength[i])
      refHeader <- CONTIGS$refHeader[i]
      CONTIG_coverage <- coverage[coverage$CHROM == refHeader, ]
      
      ################################################################################
      ##  GC normalization
      if(GCnorm == TRUE){coverage$COVERAGE <- gcNormalization(coverage$COVERAGE, fastaFile, refHeader)}
      ################################################################################
      ##  MAPPABILITY normalization
      if(MAPnorm == TRUE & file.exists(bedgraphFile)){coverage$COVERAGE <- mappabilityNormalization(coverage$COVERAGE, bedgraphFile, refHeader)}
      ################################################################################
      ##  oriC normalization
      if(ORICnorm == TRUE & !is.null(oriCposition)){coverage$COVERAGE <- oricNormalization_v2(coverage$COVERAGE, referenceLength, oriCposition)}
      ################################################################################
      ## CNV EVENTS
      peaks <- cnvOutliers(coverage, cores, peakDistanceThreshold)
      ################################################################################
      ## CNV EVENTS BORDERS DETECTION based on slope of a peak
      CNV_DF <- cnvBoundaries(peaks, coverage, step, TRIGGER)
      ################################################################################
      ## MERGE
      CNV_DF <- subset(CNV_DF, LENGTH >= 1) # DELETE SINGLE BASE EVENTS 
      ################################################################################
      # DISCORDANT READS
      CNV_DF <- cnvDiscordantReads(bamFile, CNV_DF, refHeader, referenceLength)
      ################################################################################
      
      
      ## ADD CNV copy number
      if(nrow(CNV_DF)>0){
        ## ADD CONTIG NAME
        CNV_DF$CHROM <- refHeader
        CNV_DF$COPY_NUMBER <- round(CNV_DF$COVERAGE/averageCoverage)
      } else{CNV_DF <- data.frame(ID=as.character(), CHROM=as.character(), START=integer(), END=integer(), LENGTH=integer(), COVERAGE=integer(), COPY_NUMBER=integer(),TYPE=character(),SUBTYPE=as.character(),
                                  reads_TOTAL=integer(),reads_supporting_Deletion=integer(),reads_supporting_Tandem_Direct=integer(),reads_supporting_Tandem_Indirect=integer(),
                                  reads_supporting_Interspersed_Direct=integer(),reads_supporting_Interspersed_Indirect=integer())
      }
      
      
      #reorder columns
      if(nrow(CNV_DF)>0){
        CNV_DF <- CNV_DF[,c("ID",
                            "CHROM",
                            "START",
                            "END",
                            "LENGTH",
                            "COVERAGE",
                            "COPY_NUMBER",
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
      
      CNV_DF <- as.data.frame(CNV_DF)
      MERGED_CNV_DF <- rbind(MERGED_CNV_DF, CNV_DF)
    }
    
    ## reformat and remove NAs CNVs
    MERGED_CNV_DF <- MERGED_CNV_DF[!(is.na(MERGED_CNV_DF$START)),] # remove NA CNVs
    CNV_DF <- as.data.table(MERGED_CNV_DF)
    CNV_DF <- unique(CNV_DF,by=c("CHROM","START","END","TYPE","SUBTYPE"))
    CNV_DF <- CNV_DF[order(CNV_DF$CHROM, CNV_DF$START),]
    
  }
  
  return(CNV_DF)
  ## end of function CNproScan
}
