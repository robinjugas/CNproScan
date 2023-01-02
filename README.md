# CNproScan
CNproScan is R package developed for CNV detection in bacterial genomes. It employs Generalized Extreme Studentized Deviate test for outliers to detect CNVs in read-depth data with discordant reads detection to annotate the CNVs type. 

The not-updated Matlab version is here: https://github.com/robinjugas/CNproScanMatlab

This is the latest version v1.0 For previous versions see Tags/Releases. 

## Dependencies:
Package was tested on R 4.x with several dependencies: parallel, foreach, doParallel, seqinr, Rsamtools, GenomicRanges, IRanges, data.table. 

## Installation
```
devtools::install_github("robinjugas/CNproScan")
```


## Input files:
Several input files are neccessary:
<ol>
<li>reference sequence FASTA file used in the read alignment</li>
<li>sorted and indexed BAM file from the read-aligner (BWA assumed) </li>

```
bwa index -a is reference.fasta
samtools faidx reference.fasta
bwa mem reference.fasta read1.fq read2.fq > file.sam
samtools view -b -F 4 file.sam > file.bam # mapped reads only
samtools sort -o file.bam file1.bam
samtools index file.bam
```

<li>coverage file (including zero values with -a swtich) </li>

```
samtools depth -a file.bam > file.coverage
```

<li>genome mappability file - obtained by GENMAP (https://github.com/cpockrandt/genmap) - only for the mappability normalization </li>

```
genmap index -F reference.fasta -I mapp_index
genmap map -K 30 -E 2 -I mapp_index -O mapp_genmap -t -w -bg
```

<li>origin of replication position/s - obtained from DoriC (https://origin.tubic.org/doric/browse/bacteria) - only for the oriC normalization </li>


</ol>

## Usage:
R script:
```
library("CNproScan")
# Working directory with files
setwd("workdir")
# File paths
fasta_file <- "reference.fasta"
bam_file <- "file.bam"
coverage_file <- "file.coverage"
bedgraph_file <- "mapp_genmap.bedgraph"

# For only GC normalization
DF <- CNproScanCNV(coverage_file, bam_file, fasta_file, 
                   GCnorm=TRUE, MAPnorm=FALSE, cores=4)
# Without any normalization
DF <- CNproScanCNV(coverage_file, bam_file, fasta_file, 
                   GCnorm=FALSE, MAPnorm=FALSE, cores=4)
# Both GC normalization and mappability normalization
DF <- CNproScanCNV(coverage_file, bam_file, fasta_file, 
                   GCnorm=TRUE, MAPnorm=TRUE, bedgraph_file, cores=4)

# Both GC normalization, mappability normalization and OriC normalization
DF <- CNproScanCNV(coverage_file, bam_file, fasta_file, 
                   GCnorm=TRUE, MAPnorm=TRUE,ORICnorm=TRUE, bedgraph_file,oriCposition=1, cores=4)
# or with multiple oriC positions
DF <- CNproScanCNV(coverage_file, bam_file, fasta_file, 
                   GCnorm=TRUE, MAPnorm=TRUE,ORICnorm=TRUE, bedgraph_file, oriCposition=c(10,5000), cores=4)
                   
OriC normalization is working only in single-chromosome mode!

# Write VCF file (additional function from the package)
writeVCF(DF, "fileName.vcf")
```

## Inputs description:
<ul>
<li>coverage_file = path to the .coverage file </li>
<li>bam_file = path to the .bam file </li>
<li>fasta_file = path to the .fasta file </li> 
<li>GCnorm = TRUE/FALSE whether to do GC bias normalization </li>
<li>MAPnorm = TRUE/FALSE whether to do mappability normalization </li>
<li> bedgraph_file = path to the bedgraph file outputed from genmap tool. </li>
<li>cores = number of threads for foreach %dopar%. Recommended default = 4, or more. Tested on 6,8 and 12 cores.  </li>
</ul>

## Outputs:
<ul>
<li>dataframe containing the detected CNVs </li>
<li>VCF file named sample_cnproscan.vcf in the working directory </li>
</ul>

## In development - will be added:
<ul>
</ul>

## Recent updates:
<ul>
<li> Tweaked CNV detection </li>
<li> Tweaked identification of CNV type </li>
<li>GC and mappability normalization modified and tweaked </li>
<li>new oriC normalization </li>
<li>VCF output</li>
<li>multi chromosome/contig support</li>
</ul>

## Note for multi chromosome/contig support:
BWA ignores the rest of FASTA header after the first whitespace. CNproScan expects all the headers to be the same. That means, the FASTA headers, BAM RNAME names and coverage file from samtools contain the same contig/chrosome names.   The package uses seqinr::read.fasta where  whole.header==FALSE crops header at the first whitespace. If this behaviour is issue, please post it as github issue. 

## Citation:
Robin Jugas, Karel Sedlar, Martin Vitek, Marketa Nykrynova, Vojtech Barton, Matej Bezdicek, Martina Lengerova, Helena Skutkova,
CNproScan: Hybrid CNV detection for bacterial genomes,
Genomics, Volume 113, Issue 5, 2021, Pages 3103-3111, ISSN 0888-7543,
https://doi.org/10.1016/j.ygeno.2021.06.040.
(https://www.sciencedirect.com/science/article/pii/S0888754321002779)
