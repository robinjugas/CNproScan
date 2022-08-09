#' MAPPABILITY normalization
#' This function normalizes the coverage values by mappability of genome regions. Requires genmap tool bedgraph input.   
#' 
#' @param coverageVECTOR dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param bedgraphFile Path to the bedgraph file outputed from genmap tool. 
#' @param refHeader name of chromosome or conrig
#' @return A vector of normalized coverage values
#' @export
#' @noRd
#' 
mappabilityNormalization <- function(coverageVECTOR, bedgraphFile, refHeader){
  ## 
  referenceLength <- length(coverageVECTOR)
  coverageVECTOR_new <- rep(0, times=referenceLength)
  med <- median(coverageVECTOR)
  
  ##  READ BEDGRAPH FILE
  mappability <- read.csv(bedgraphFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(mappability) != 4){
    stop("`bedgraphFile` must have 4 columns as outputed by 'genmap'. ")
  }
  header <- c('ID', 'START','STOP','MAPPABILITY')
  colnames(mappability) <- header
  mappability <- mappability[mappability$ID == refHeader, ]
  mappability <- mappability[, -1] #keeps only POS and coverageVECTOR
  # BEDgraph is zero based
  mappability[, 'START'] <-  mappability[, 'START']+1
  
  ## Calculate READ count in BEDGRAPH regions
  RC_count <- rep(0, times=nrow(mappability))
  for (i in 1:nrow(mappability)){
    RC_count[i] <- mean( coverageVECTOR[mappability[i,"START"]:mappability[i,"STOP"]])
  }
  mappability[, 'READDEPTH'] <- RC_count
  mappability[, 'MAPPABILITY'] <- round(mappability[, 'MAPPABILITY']*100)
  # median MAPPscore
  medReadCount <- median(RC_count) #median RC for all bins
  
  ## match unique mappability scores with median of read-depths means in those regions with same mappability
  #MAPPmat = 1.col - unique mappability scores 2.col - median coverageVECTOR of all regions with given mappability score
  uniqueMAPPscore <- sort(unique(mappability[, 'MAPPABILITY']))
  MAPPmat <- data.frame(MAPPABILITY=integer(), MED_COV=integer())
  MAPPmat[1:length(uniqueMAPPscore), "MAPPABILITY"] <- uniqueMAPPscore
  for (i in 1:nrow(MAPPmat)){
    value <- MAPPmat[i, "MAPPABILITY"]
    tempindices <- which(mappability[, "MAPPABILITY"] == value)
    MAPPmat[i, "MED_COV"] <- median(mappability[tempindices, "READDEPTH"])
  }
  
  ## find corresponding median read-depth of given Mappability score for every region
  for  (i in 1:nrow(mappability)){
    findices <- which(MAPPmat[, "MAPPABILITY"] == mappability[i, "MAPPABILITY"])
    mappability[i, "MEDIANRC"] <- MAPPmat[findices, "MED_COV"]
  }
  
  ## Normalization
  for  (i in 1:nrow(mappability)){
    start <- mappability[i, "START"]
    stop <- mappability[i, "STOP"]
    coeficient <- medReadCount/mappability[i, "MEDIANRC"]
    coverageVECTOR_new[start:stop] <- coverageVECTOR[start:stop]*coeficient
    
  }
  # diff <- abs(coverageVECTOR_new-coverageVECTOR)
  # plot(diff)
  return(coverageVECTOR_new)
}

