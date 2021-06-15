#' MAPPABILITY normalization
#' This function normalizes the coverage values and limits the GC content bias induced by PCR.  
#' Inputs: coverage dataframe (1st col POS, 2nd col CCOVERAGE), path to bedgraph file
#' Output: vector of modified COVERAGE values
#' 
#' @param infile Path to the input file
#' @return A vector of new coverage values
#' @export
#' 
MAPPABILITYnormalization <- function(coverage, bedgraphFile){
  ## 
  referenceLength <- length(coverage)
  coverage_new <- rep(0, times=referenceLength)
  med <- median(coverage)
  
  ##  READ BEDGRAPH FILE
  mappability <- read.csv(bedgraphFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  header <- c('ID', 'START','STOP','MAPPABILITY')
  colnames(mappability) <- header
  mappability <- mappability[,-1] #keeps only POS and COVERAGE
  # BEDgraph is zero based
  mappability[,'START'] <-  mappability[,'START'] +1
  
  ## Calculate READ count in BEDGRAPH regions
  RC_count <- rep(0, times=nrow(mappability))
  for (i in 1:nrow(mappability)){
    RC_count[i] <- mean( coverage[mappability[i,"START"]:mappability[i,"STOP"]])
  }
  mappability[,'READDEPTH'] <- RC_count
  mappability[,'MAPPABILITY'] <- round(mappability[,'MAPPABILITY']*100)
  # median MAPPscore
  medReadCount=median(RC_count) #median RC for all bins
  
  ## match unique mappability scores with median of read-depths means in those regions with same mappability
  #MAPPmat = 1.col - unique mappability scores 2.col - median coverage of all regions with given mappability score
  uniqueMAPPscore <- sort(unique(mappability[,'MAPPABILITY']))
  MAPPmat <- data.frame(MAPPABILITY=integer(),MED_COV=integer())
  MAPPmat[1:length(uniqueMAPPscore),"MAPPABILITY"] <- uniqueMAPPscore
  for (i in 1:nrow(MAPPmat)){
    value <- MAPPmat[i,"MAPPABILITY"]
    tempindices <- which(mappability[,"MAPPABILITY"]==value)
    MAPPmat[i,"MED_COV"] <- median(mappability[tempindices,"READDEPTH"])
  }
  
  ## find corresponding median read-depth of given Mappability score for every region
  for  (i in 1:nrow(mappability)){
    findices <- which(MAPPmat[,"MAPPABILITY"]==mappability[i,"MAPPABILITY"])
    mappability[i,"MEDIANRC"] <- MAPPmat[findices,"MED_COV"]
  }
  
  ## Normalization
  for  (i in 1:nrow(mappability)){
    start <- mappability[i,"START"]
    stop <- mappability[i,"STOP"]
    coeficient <- medReadCount/mappability[i,"MEDIANRC"]
    coverage_new[start:stop] <- coverage[start:stop]*coeficient
    
  }
  # diff <- abs(coverage_new-coverage)
  # plot(diff)
  return(coverage_new)
}

