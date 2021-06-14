# CNproScan
CNproScan is R package developed for CNV detection in bacterial genomes. It employs Generalized Extreme Studentized Deviate test for outliers to detect CNVs in read-depth data with discordant reads detection to annotate the CNVs. 
The package is still in development. Namely, the normalization part of package is still in development. 

## Dependencies:
Package was tested on R 4.1.0 with several dependencies: parallel, foreach, doParallel, seqinr, Rsamtools, GenomicRanges, IRanges. 

## Installation
```
devtools::install_github("robinjugas/CNproScan")
```


## Input files:
Several input files are neccessary:
<ol>
<li>reference sequence FASTA file used in the read alignment</li>
<li>sorted and indexed BAM file from the read-aligner </li>
```
bwa index -a is reference.fasta
samtools faidx reference.fasta
bwa mem reference.fasta read1.fq read2.fq > file.sam
samtools view -b -F 4 file.sam > file.bam # mapped reads only
samtools sort -o file.bam file1.bam
samtools index file.bam
```

<li>coverage file (including zero values) </li>
```
samtools depth -a file.bam > file.coverage
```

<li>genome mappability file - obtained by GENMAP (https://github.com/cpockrandt/genmap) </li>
```
genmap index -F reference.fasta -I mapp_index
genmap map -K 30 -E 2 -I mapp_index -O mapp_genmap -t -w -bg
```

<li>oriC position - use DoriC database (http://tubic.org/doric/public/index.php/search)</li>
</ol>

## Usage:
R script:
```
CNproScanCNV(coverage_file, bam_file, fasta_file,number of threads)
coverage_file = path to the .coverage file
bam_file = path to the .bam file
fasta_file = 
number of threads = 


```

## In development - will be added:
<ul>
<li>normalization of read-depth</li>
<li>VCF output</li>
</ul>