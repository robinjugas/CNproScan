#' CNproScan - CNV detection for bacteria genomes
#' This function detects and annotate CNVs in bacterial genomes. It uses GESD outliers detection 
#' and detection of discordant read-pairs for annotation. 
#' 
#' @importFrom seqinr getLength read.fasta
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
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

CNproScanCNV <- function(coverageFile,bamFile,fastaFile,GCnorm=TRUE,MAPnorm=FALSE,bedgraphFile=NULL,cores=2){
  
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
    # coverageORI   <- oricNormalization(coverage, oriC_position=517)
    ################################################################################
    ## DETECTING DELETIONS AS ZERO COVERAGE
    TEMP <- cnvDeletions(coverage, peakDistanceThreshold)
    coverage <- TEMP[[1]]
    DEL_DF <- TEMP[[2]]
    ################################################################################
    ## CNV EVENTS
    peaks <- cnvOutliers(coverage, cores, peakDistanceThreshold)
    ################################################################################
    ## CNV DUPLICATIONS EVENTS BORDERS DETECTION based on slope of a peak
    DUP_DF <- cnvBoundaries(peaks, coverage, averageCoverage, step, TRIGGER)
    ################################################################################
    ## MERGE
    CNV_DF <- rbind(DEL_DF, DUP_DF) # MERGE DELETIONS AND DUPLICATIONS
    CNV_DF <- subset(CNV_DF, LENGTH >= 1) # DELETE SINGLE BASE EVENTS 
    CNV_DF[CNV_DF$COVERAGE < averageCoverage, "TYPE"] <- "DEL" #CHECK FOR CNV < averageCoverage - those are likely deletions
    gc()
    ################################################################################
    # DISCORDANT READS
    CNV_DF <- cnvDiscordantReads(bamFile, CNV_DF, refHeader, referenceLength)
    ################################################################################
    ## ADD CONTIG NAME
    CNV_DF$CHROM <- refHeader
    
    col_order <- c("ID", "CHROM", colnames(CNV_DF)[3:ncol(CNV_DF)-1])
    CNV_DF <- CNV_DF[, col_order]
    
    ## write VCF
    # sampleName <- strsplit(basename(coverageFile), "[.]")[[1]][1]
    # writeVCF(CNV_DF, sampleName)
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
      if(GCnorm == TRUE){CONTIG_coverage$COVERAGE <- gcNormalization(CONTIG_coverage$COVERAGE, fastaFile, refHeader)}
      ################################################################################
      ##  MAPPABILITY normalization
      if(MAPnorm == TRUE & file.exists(bedgraphFile)){CONTIG_coverage$COVERAGE <- mappabilityNormalization(CONTIG_coverage$COVERAGE, bedgraphFile, refHeader)}
      ################################################################################
      ##  oriC normalization
      ## coverageORI   <- oricNormalization(CONTIG_coverage, oriC_position=517)
      ################################################################################
      ## DETECTING DELETIONS AS ZERO COVERAGE
      TEMP <- cnvDeletions(CONTIG_coverage, peakDistanceThreshold)
      CONTIG_coverage <- TEMP[[1]]
      DEL_DF <- TEMP[[2]]
      ################################################################################
      ## CNV EVENTS
      peaks <- cnvOutliers(CONTIG_coverage, cores, peakDistanceThreshold)
      ################################################################################
      ## CNV DUPLICATIONS EVENTS BORDERS DETECTION based on slope of a peak
      DUP_DF <- cnvBoundaries(peaks, CONTIG_coverage, averageCoverage, step, TRIGGER)
      ################################################################################
      ## MERGE
      CNV_DF <- rbind(DEL_DF, DUP_DF) # MERGE DELETIONS AND DUPLICATIONS
      CNV_DF <- subset(CNV_DF, LENGTH>=1) # DELETE SINGLE BASE EVENTS 
      CNV_DF[CNV_DF$COVERAGE<averageCoverage, "TYPE"] <- "DEL" #CHECK FOR CNV < averageCoverage - those are likely deletions
      gc()
      ################################################################################
      ## DISCORDANT READS
      CNV_DF <- cnvDiscordantReads(bamFile, CNV_DF, refHeader, referenceLength)
      ## ADD CONTIG NAME
      CNV_DF$CHROM <- refHeader
      
      col_order <- c("ID", "CHROM", colnames(CNV_DF)[3:ncol(CNV_DF)-1])
      CNV_DF <- CNV_DF[, col_order]
      
      MERGED_CNV_DF <- rbind(MERGED_CNV_DF, CNV_DF)
    }
    
    ## write VCF
    CNV_DF <- MERGED_CNV_DF
    CNV_DF <- CNV_DF[order(CNV_DF$CHROM, CNV_DF$START),]
    # sampleName <- strsplit(basename(coverageFile), "[.]")[[1]][1]
    # writeVCF(CNV_DF, sampleName)
  }
  
  
  return(CNV_DF)
  ## end of function CNproScan
}
