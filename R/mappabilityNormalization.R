#' WINDOWS normalization v2
#' This function normalizes the coverage values by mappability of genome regions. Requires genmap tool bedgraph input.   
#' 
#' @param coverageVECTOR dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param bedgraphFile Path to the bedgraph file outputed from genmap tool. 
#' @param refHeader name of chromosome or contig
#' @return A vector of normalized coverage values
#' @export
#' @noRd
#' 
mappabilityNormalization <- function(coverageVECTOR, bedgraphFile, refHeader){
  ## 
  coverageVECTOR_new <- coverageVECTOR

  ##  READ BEDGRAPH FILE
  WINDOWS <- read.csv(bedgraphFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(WINDOWS) != 4){
    stop("`bedgraphFile` must have 4 columns as outputed by 'genmap'. ")
  }
  header <- c('ID', 'START','STOP','MAPPABILITY')
  colnames(WINDOWS) <- header
  WINDOWS <- WINDOWS[WINDOWS$ID == refHeader, ]
  WINDOWS <- WINDOWS[, -1] #keeps only POS and coverageVECTOR
  # BEDgraph is zero based
  WINDOWS[, 'START'] <-  WINDOWS[, 'START']+1
  
  ## Calculate READ count in BEDGRAPH regions
  RC_count <- rep(0, times=nrow(WINDOWS))
  for (i in 1:nrow(WINDOWS)){
    RC_count[i] <- mean( coverageVECTOR[WINDOWS[i,"START"]:WINDOWS[i,"STOP"]])
  }
  WINDOWS[, 'READDEPTH'] <- round(RC_count)
  WINDOWS[, 'MAPPABILITY'] <- round(WINDOWS[, 'MAPPABILITY']*100)
  # median MAPPscore
  medReadCount <- median(RC_count) #median RC for all bins
  
  ## match unique WINDOWS scores with median of read-depths means in those regions with same WINDOWS
  ## MAPP_tab = 1.col - unique WINDOWS scores 2.col - median coverageVECTOR of all regions with given WINDOWS score
  uniqueMAPPscore <- sort(unique(WINDOWS[, 'MAPPABILITY']))
  MAPP_tab <- data.frame(MAPPABILITY=integer(), MED_COV=integer())
  MAPP_tab[1:length(uniqueMAPPscore), "MAPPABILITY"] <- uniqueMAPPscore
  for (i in 1:nrow(MAPP_tab)){
    tempindices <- which(WINDOWS[, "MAPPABILITY"] == MAPP_tab[i, "MAPPABILITY"])
    MAPP_tab[i, "MED_COV"] <- median(na.omit(WINDOWS$READDEPTH[tempindices]))
  }
  
  ## find corresponding median read-depth of given WINDOWS score for every region
  for  (i in 1:nrow(WINDOWS)){
    findices <- which(MAPP_tab[, "MAPPABILITY"] == WINDOWS[i, "MAPPABILITY"])
    WINDOWS[i, "MEDIANRC"] <- MAPP_tab[findices, "MED_COV"]
  }
  
  ## Mappability Normalization
  for  (i in 1:nrow(WINDOWS)){
    start <- WINDOWS[i, "START"]
    stop <- WINDOWS[i, "STOP"]
    
    if(WINDOWS[i, "MEDIANRC"]>1){
      coeficient <- medReadCount/WINDOWS[i, "MEDIANRC"]
      coverageVECTOR_new[start:stop] <- coverageVECTOR_new[start:stop]*coeficient
    } # no normalization for extreme values
  }

  return(coverageVECTOR_new)
  
}

