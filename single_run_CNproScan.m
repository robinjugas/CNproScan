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

