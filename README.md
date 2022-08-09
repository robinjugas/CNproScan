# CNproScan
CNproScan is R package developed for CNV detection in bacterial genomes. It employs Generalized Extreme Studentized Deviate test for outliers to detect CNVs in read-depth data with discordant reads detection to annotate the CNVs type. 

The not-updated Matlab version is here: https://github.com/robinjugas/CNproScanMatlab


## Dependencies:
Package was tested on R 4.x with several dependencies: parallel, foreach, doParallel, seqinr, Rsamtools, GenomicRanges, IRanges. 

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

</ol>

## Usage:
R script:
```
CNproScanCNV(coverage_file, bam_file, fasta_file,number of threads)
CNproScanCNV(coverageFile,bamFile,fastaFile,GCnorm=TRUE,MAPnorm=FALSE,cores=12)
```
Inputs:
coverageFile = path to the .coverage file
bamFile = path to the .bam file
fastaFile = path to the .fasta file
GCnorm = TRUE/FALSE whether to do GC bias normalization
MAPnorm = TRUE/FALSE whether to do mappability normalization
cores = number of threads for foreach %dopar%. Recommended default = 4, or more. Tested on 6,8 and 12 cores.  

Outputs:
dataframe containing the detected CNVs
VCF file named sample_cnproscan.vcf in the working directory. 


## In development - will be added:
<ul>
<li>oriC normalization </li>
</ul>

## Recent updates:
<ul>
<li>GC and mappability normalization </li>
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
