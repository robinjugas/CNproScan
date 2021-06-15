#' GC normalization
#' This function normalizes the coverage values and limits the GC content bias induced by PCR.  
#' Inputs: coverage dataframe (1st col POS, 2nd col CCOVERAGE), path to fasta file
#' Output:  vector of modified COVERAGE values
#' 
#' @importFrom seqinr getLength read.fasta GC
#' @param infile Path to the input file
#' @return A vector of new coverage values
#' @export
#' 
GCnormalization <- function(coverage, fastaFile){
  ## 
  windowLength <- 100
  referenceLength <- length(coverage)
  GC_count <- rep(0, times=referenceLength)
  coverage_new <- rep(0, times=referenceLength)
  med <- median(coverage)
  
  ##  READ FASTA FILE
  reference <- seqinr::read.fasta(fastaFile,seqonly = TRUE)
  reference <- reference[[1]]
  reference <- strsplit(reference, "")[[1]]
  ## Calculate GC content in window
  for (i in seq(from = 1, to = referenceLength, by = windowLength/2)){
    if (i <= (referenceLength-windowLength)){
      GC_count[i:(i+windowLength-1)] <- seqinr::GC(reference[i:(i+windowLength-1)],NA.GC=NA)
    }
    else{
      GC_count[i:length(GC_count)] <- seqinr::GC(reference[i:referenceLength],NA.GC=NA)
    }
  }
  GC_count <- GC_count*100 #make a percentage
  
  ## GC dataframe GC_tab - 1.col - GC value in % 2.col - median coverage of all windows with given GC % value
  GC_tab <- data.frame(GC=integer(),MED_COV=integer())
  GCsek <- seq(from=min(GC_count), to=max(GC_count))
  GC_tab[1:length(GCsek),"GC"] <- GCsek
  GC_tab[,"MED_COV"] <- rep(0, times=length(GCsek))
  for(i in 1:nrow(GC_tab)){
    position <- which(GC_count==GC_tab[i, "GC"])
    pom <- median(coverage[position])
    GC_tab[i,"MED_COV"]  <- pom
  }
  
  ## New normalized coverage
  for(i in 1:referenceLength){
    GC_act <- GC_count[i]
    if(GC_act==0 | GC_tab[GC_act,"GC"]==0 | is.na(GC_tab[GC_act,"MED_COV"]) | coverage[i]<5 ){
      coverage_new[i] <- coverage[i] # zustava puvodni pokryti
    }
    else
    {
      med2 <- GC_tab[GC_act,"MED_COV"] # vyhledani medianu pozic pokryti s danym GC poctem ve vytvorene tabulce
      coverage_new[i] <- coverage[i]*(med/med2) 
    }
  }
  # diff <- abs(coverage_new-coverage)
  # plot(diff)
  return(coverage_new)
}

