# CNproScan
CNproScan is a set of MATLAB functions to detect CNV in bacteria genomes. 

## Installation
Clone the whole repository or the folder with MATLAB files (CNV_detekce_final_v6_CNV2) and modify the running script. 

## Input files:
Several input files are neccessary:

1. reference FASTA file used in the read alignment

2. BAM file from read-aligner
```
bwa index -a is reference.fasta
samtools faidx reference.fasta
bwa mem reference.fasta read1.fq read2.fq > file.sam

samtools view -b -F 4 file.sam > file.bam # mapped reads only
samtools sort -o file.bam file1.bam
samtools index file.bam
```
3. coverage file
```
samtools depth -a file.bam > file.coverage
```
4. genome mappability file - obtained by GENMAP (https://github.com/cpockrandt/genmap)
```
genmap index -F reference.fasta -I mapp_index
genmap map -K 30 -E 2 -I mapp_index -O mapp_genmap -t -w -bg
```

5. oriC position - use DoriC database (http://tubic.org/doric/public/index.php/search)

## Usage:
MATLAB script:
```
clc, clear all

addpath('CNV_detekce_final_v6_CNV2'); %folder with CNproScan functions

%% INITIALIZATION
coverage_file = 'CNVseq.coverage';
reference_file = 'FN433596.fasta';
genmap_file = 'FN433596.bedgraph';
bamfile = 'CNVseq.sorted.bam';
oriC_position=517;
step=100;

%% COMPUTING

% Loading read-depth files
coverage=CNV1_import_coverage(coverage_file, reference_file, genmap_file, oriC_position);
coverageSignal=coverage(:,2)'; %take only read-depth values

% Peaks detection
[peaksPolished,indicationPeaks]=CNV2_peaks_detection(coverageSignal);

% Read-pairs distance detection
[distanceSignal,insertSize]=CNV3_pair_reads_distance(bamfile,step);

% Read-pairs distance thresholding
[indicationHigher,indicationLower]=CNV4_pair_reads_distance_thresholding(distanceSignal,insertSize);

% Final output
[CNVtable,CNVseq] = CNV5_detection_output(peaksPolished,indicationPeaks,indicationHigher,indicationLower,coverageSignal,reference_file);

% Writing into excel spreasheet
writecell(CNVtable,['CNV_detection_v6.xls'])
```
