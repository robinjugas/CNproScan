#' gcNormalization v2
#' This function normalizes the coverage values and limits the GC content bias induced by PCR.  
#' Inputs: coverage dataframe (1st col POS, 2nd col CCOVERAGE), path to fasta file
#' Output:  vector of modified COVERAGE values
#' 
#' @importFrom seqinr getLength read.fasta GC
#' @param coverageVECTOR dataframe of coverage, 3 columns CHROM POS COVERAGE
#' @param fastaFile path to the FASTA file
#' @param refHeader name of chromosome or conrig
#' @return A vector of new coverage values
#' @export
#' @noRd
#' 
gcNormalization <- function(coverageVECTOR, fastaFile, refHeader){
  ## 
  windowLength <- 100
  referenceLength <- length(coverageVECTOR)
  coverageVECTOR_new <- coverageVECTOR
  
  ##  READ FASTA FILE
  reference <- seqinr::read.fasta(fastaFile, seqonly = FALSE, whole.header=FALSE)
  reference <- reference[[which(names(reference)==refHeader)]]

  referenceLength_contig <- length(reference)

  ## DATA.FRAME with windows
  WINDOWS <- data.frame(START=integer(), STOP=integer(),GC=integer(),READDEPTH=integer())
  WINDOWS[1:length(seq(from=1, to=referenceLength_contig, by=windowLength)),"START"] <- seq(from=1, to=referenceLength_contig, by=windowLength)
  WINDOWS$STOP <- WINDOWS$START+windowLength-1
  WINDOWS <- WINDOWS[-nrow(WINDOWS)] #delete last window
  WINDOWS$STOP[nrow(WINDOWS)] <- referenceLength_contig #extend the last windows to the end
  
  ## Calculate GC content  and READ count in window
  RC_count <- rep(0, times=nrow(WINDOWS))
  GC_perc <- rep(0, times=nrow(WINDOWS))
  for (i in 1:nrow(WINDOWS)){
    RC_count[i] <- mean(coverageVECTOR[WINDOWS[i,"START"]:WINDOWS[i,"STOP"]])
    GC_perc[i] <- seqinr::GC(reference[WINDOWS[i,"START"]:WINDOWS[i,"STOP"]], NA.GC=0)
    
  }
  WINDOWS[, "READDEPTH"] <- round(RC_count)
  WINDOWS[, "GC"] <- round(GC_perc*100)
  
  ## median read depth of windows
  medReadCount <- median(RC_count) #median RC for all bins
  
  ## GC dataframe GC_tab - 1.col - GC value in % 2.col - median coverageVECTOR of all windows with given GC % value
  GC_tab <- data.frame(GC=integer(), MED_COV=integer())
  GCsek <- seq(from=0, to=100)
  GC_tab[1:length(GCsek), "GC"] <- GCsek
  GC_tab[, "MED_COV"] <- rep(0, times=length(GCsek))
  for(i in 1:nrow(GC_tab)){
    position <- which(WINDOWS$GC==GC_tab[i, "GC"])
    GC_tab[i, "MED_COV"]  <- median(na.omit(WINDOWS$READDEPTH[position]))
  }
  
  ## find corresponding median read-depth of given GC percentage score for every region
  for  (i in 1:nrow(WINDOWS)){
    idx <- which(GC_tab[, "GC"] == WINDOWS[i, "GC"])
    WINDOWS[i, "MEDIANRC"] <- GC_tab[idx, "MED_COV"]
  }
  
  
  ## GC Normalization
  for  (i in 1:nrow(WINDOWS)){
    start <- WINDOWS[i, "START"]
    stop <- WINDOWS[i, "STOP"]
    
    if(WINDOWS[i, "MEDIANRC"]>1){
      coeficient <- medReadCount/WINDOWS[i, "MEDIANRC"]
      coverageVECTOR_new[start:stop] <- coverageVECTOR_new[start:stop]*coeficient
    }# no normalization for extreme values
    
  }
  
  return(coverageVECTOR_new)
  
}

