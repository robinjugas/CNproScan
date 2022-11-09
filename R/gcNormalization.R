#' gcNormalization
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
  GC_count <- rep(0, times=referenceLength)
  coverageVECTOR_new <- rep(0, times=referenceLength)
  med <- median(coverageVECTOR)
  
  ##  READ FASTA FILE
  reference <- seqinr::read.fasta(fastaFile, seqonly = FALSE, whole.header=FALSE)
  # refHeader %in% names(reference) 
  # which(names(reference)==refHeader)
  reference <- reference[[which(names(reference)==refHeader)]]
  
  
  ## Calculate GC content in window
  for (i in seq(from=1, to=referenceLength, by=windowLength/2)){
    if (i <= (referenceLength-windowLength)){
      GC_count[i:(i+windowLength-1)] <- seqinr::GC(reference[i:(i+windowLength-1)], NA.GC=NA)
    }
    else{
      GC_count[i:length(GC_count)] <- seqinr::GC(reference[i:referenceLength], NA.GC=NA)
    }
  }
  GC_count <- GC_count*100 #make a percentage
  
  ## GC dataframe GC_tab - 1.col - GC value in % 2.col - median coverageVECTOR of all windows with given GC % value
  GC_tab <- data.frame(GC=integer(), MED_COV=integer())
  GCsek <- seq(from=0, to=100)
  GC_tab[1:length(GCsek), "GC"] <- GCsek
  GC_tab[, "MED_COV"] <- rep(0, times=length(GCsek))
  for(i in 1:nrow(GC_tab)){
    position <- which(GC_count==GC_tab[i, "GC"])
    pom <- median(coverageVECTOR[position])
    GC_tab[i, "MED_COV"]  <- pom
  }
  
  ## New normalized coverage
  for(i in 1:referenceLength){
    GC_act <- as.integer(GC_count[i])
    if( GC_act == 0 | GC_tab[which(GC_tab$GC==GC_act), "GC"] == 0 | is.na(GC_tab[which(GC_tab$GC==GC_act), "MED_COV"]) | coverageVECTOR[i] < 1  ){ 
      coverageVECTOR_new[i] <- coverageVECTOR[i] # zustava puvodni pokryti
    }else{
      med2 <- GC_tab[which(GC_tab$GC==GC_act), "MED_COV"] # vyhledani medianu pozic pokryti s danym GC poctem ve vytvorene tabulce
      coverageVECTOR_new[i] <- coverageVECTOR[i]*(med/med2) 
    }
  }
  
  # ##COMPARISON
  # diff <- abs(coverageVECTOR_new-coverageVECTOR)
  # plot(diff)
  # #original dependence
  # plot(x=GC_tab$GC,y=GC_tab$MED_COV)
  # ## GC dataframe GC_tab - 1.col - GC value in % 2.col - median coverageVECTOR of all windows with given GC % value
  # GC_tab <- data.frame(GC=integer(), MED_COV=integer())
  # GCsek <- seq(from=0, to=100)
  # GC_tab[1:length(GCsek), "GC"] <- GCsek
  # GC_tab[, "MED_COV"] <- rep(0, times=length(GCsek))
  # for(i in 1:nrow(GC_tab)){
  #   position <- which(GC_count==GC_tab[i, "GC"])
  #   pom <- median(coverageVECTOR_new[position])
  #   GC_tab[i, "MED_COV"]  <- pom
  # }
  # #new dependence
  # plot(x=GC_tab$GC,y=GC_tab$MED_COV)
  # plot(coverageVECTOR)
  # plot(coverageVECTOR_new)
  
  ##############################################
  
  return(coverageVECTOR_new)
}

