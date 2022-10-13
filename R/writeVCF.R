#' Write VCF file from CNproScan dataframe
#' Function to reformat CNV dataframe to VCF formated dataframe to write down
#' Inputs: detected CNV dataframe 
#' Output: VCF file
#' 
#' @param CNV_DF CNV dataframe
#' @param fileName name of the VCF file
#' @return VCF dataframe ready to write to disk
#' @export
#' 
writeVCF <- function(CNV_DF,fileName){

  VCF_tab <- data.frame(CHROM=character(), POS=character(), ID=character(), REF=character(), ALT=character(), QUAL=character(),
                        FILTER=character(), INFO=character(), FORMAT=character(), SAMPLE=character())
  
  # write VCFheader
  vcf_header <- read.delim(system.file('extdata', "CNproscan_vcf_header.txt", package='CNproScan'))
  
  file.create(fileName, overwrite=TRUE)
  
  for(i in 1:nrow(vcf_header)){ 
    write(vcf_header[[1]][i],  file=fileName, append=TRUE)
  }
  
  header <- paste("#CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", sep="\t")
  write(header, file=fileName, append=TRUE)
  
  # CNV_DF <- unique(CNV_DF,by=c("CHROM","START","END","TYPE"))
  
  for(i in 1:nrow(CNV_DF)){
    VCF_tab[i, "CHROM"] <- CNV_DF[i,"CHROM"] # chromosome
    VCF_tab[i, "POS"] <- CNV_DF[i,"START"] #start
    VCF_tab[i, "ID"] <- paste0(i,"_cnproscan") # some ID, number of CNV
    VCF_tab[i, "REF"] <- "N" #lumpy inspired
    VCF_tab[i, "ALT"] <- paste0("<",CNV_DF[i,"TYPE"],">") #lumpy inspired
    VCF_tab[i, "QUAL"] <- "." #lumpy inspired
    VCF_tab[i, "FILTER"] <- "PASS" #cnvnator inspired
    
    VCF_tab[i, "INFO"] <- paste(paste0("END=",CNV_DF[i,"END"]),paste0("SVTYPE=",CNV_DF[i,"TYPE"]),paste0("SVLEN=",CNV_DF[i,"LENGTH"]),paste0("SVSUBTYPE=",CNV_DF[i,"SUBTYPE"]),sep=";")
    VCF_tab[i, "FORMAT"] <- "."
    VCF_tab[i, "SAMPLE"] <- "."
    
    temp <- do.call(paste, c(VCF_tab[c(1:ncol(VCF_tab))], sep="\t"))
    write(temp, file=fileName, append=TRUE)
  }

}

##fileformat=VCFv4.1
##fileDate=20220222
##reference=KP_ref.fasta
##source=CNVnator
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">
##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">
##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">
##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">
##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">
##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">
##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">
##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	results/cnvnator/S02_L001/S02_L001
# gi|238892256|ref|NC_012731.1|	17801	S02_L001_CNVnator_del_1	N	<DEL>	.	PASS	END=18000;SVTYPE=DEL;SVLEN=-200;IMPRECISE;natorRD=0.0486348;natorP1=165101;natorP2=5.85755e-06;natorP3=1;natorP4=1;natorQ0=0.0555556	GT:CN	1/1:0
# 
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S02_L001
# gi|238892256|ref|NC_012731.1|	122141	1	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+
#   
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S02_L001
# gi|238892256|ref|NC_012731.1|	122141	1	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+:28;SVLEN=92084;END=214225;CIPOS=-130,7;CIEND=-10,105;CIPOS95=-15,1;CIEND95=-2,12;IMPRECISE;SU=28;PE=28;SR=0	GT:SU:PE:SR	./.:28:28:0
