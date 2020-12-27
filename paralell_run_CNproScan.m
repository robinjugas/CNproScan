clear all, clc

addpath('F:\SHARED\CNV_detekce_final_v5'); %folder with CNproScan functions

%% INITIALIZATION
coverageFiles = dir('*.coverage'); % list coverage files in current directory
numfiles = length(coverageFiles); % number of files
bamFiles = dir('*.bam'); %  list bam files in current directory

reference_file = 'FN433596.fasta';
genmap_file = 'FN433596.bedgraph';
oriC_position=517;
step=100;

%% COMPUTING
j=1; % only one file, delete if multiple files tested

parpool('local',2)

parfor i = 1:2
    if i == 1
        % Loading read-depth files
        coverage=CNV1_import_coverage(coverageFiles(j,1).name, reference_file, genmap_file, oriC_position);
        coverageSignal=coverage(:,2)'; %take only read-depth values, 2nd column

        % Peaks detection
        [peaksPolished,indicationPeaks]=CNV2_peaks_detection(coverageSignal);
        
        parsave('data1.mat',coverageSignal,peaksPolished,indicationPeaks)
        
        disp('PEAKS DETECTED')
    else
        % Read-pairs distance detection
        [distanceSignal,insertSize]=CNV3_pair_reads_distance(bamFiles(j,1).name,step);
        % Read-pairs distance thresholding
        [indicationHigher,indicationLower]=CNV4_pair_reads_distance_thresholding(distanceSignal,insertSize);
        
        parsave('data2.mat',distanceSignal,indicationHigher,indicationLower)
        
        disp('READ-PAIRS SCANNED')
        
    end
end


load('data1.mat')
coverageSignal=x;
peaksPolished=y;
indicationPeaks=z;

load('data2.mat')
distanceSignal=x;
indicationHigher=y;
indicationLower=z;


if ~isempty(peaksPolished)
    % Final output
    [CNVtable,CNVseq] = CNV5_detection_output(peaksPolished,indicationPeaks,indicationHigher,...
        indicationLower,coverageSignal,FASTAfile);
    writecell(CNVtable,['CNproScan_detection_output.xls'])
else
    
    
end
delete(gcp)

% save variables in parpool function
function parsave(fname, x,y,z)
save(fname, 'x', 'y','z')
end


