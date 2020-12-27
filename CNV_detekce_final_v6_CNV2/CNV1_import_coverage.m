function [coverage_final] = CNV1_import_coverage(coverage_file,reference_file,genmap_file, oriC_position)
%% IMPORTING COVERAGE FROM bam file
% coverage_file coverage file name, obtained by samtools coverage
% samtools depth -a $file > $file.coverage
% reference_file FASTA reference file name
% genmap_file BEDGRAPH file,

%% FASTA reference
[~,ref_seq]=fastaread(reference_file); %length of reference used to BWA
ref_length=length(ref_seq);

%% COVERAGE files, extracting data by 10k rows
fileID = fopen(coverage_file);
formatSpec = '%*s %d %d';
coverage_raw = zeros(0, 2);
k = 0;
while ~feof(fileID)
    k = k+1;
    temp = textscan(fileID,formatSpec,10000);
    coverage_raw = [coverage_raw; [temp{:}]];
end
fclose(fileID);
disp('CNV1 Loading coverage DONE')

%% COMPUTING GC content in sliding overlapping window based on consensual sequence
W=70; %length of windows to be EVEN
L=ref_length;
GC=zeros(1,L);

coverage_GC=coverage_raw;

%START (windows half)
if sum(coverage_GC(1:W/2,2))<2*W %prázdné okno malá coverage
    GC(1,1:W/2)=0;
else
    bases=basecount(ref_seq(1:W/2));
    GC(1,1:W/2)=round((bases.C+bases.G)/(W/2)*100);
end

%STOP (windows half)
if sum(coverage_GC(L-W/2:L,2))<2*W %prázdné okno malá coverage
    GC(L-W/2:L)=0;
else
    bases=basecount(ref_seq(L-W/2+1:L));
    GC(L-W/2:L)=round((bases.C+bases.G)/(W/2)*100);
end

%REST of coverage
for i=W/2+1:L-W/2
    if coverage_GC(i,2)<5 % IF coverage at position <5, I skip counting GC there
        GC(i)=0;
    else
        bases=basecount(ref_seq(i-W/2:i+W/2-1));
        GC(i)=round(((bases.C+bases.G)/W)*100);
    end
end
coverage_GC(:,3)=GC'; % GC content 3rd column in cell matrix
disp('CNV1 GC content computing DONE')



%% NORMALIZATION OF GC CONTENT

% median RCs of given GCcontent
coverage_GC_normalized=coverage_GC;
GCmat=zeros(101,2);
GCmat(:,1)=[0:100]'; %1sl GC%, 2sl median Read Counts with corresponding GC content value
temp=coverage_GC_normalized(:,2:3); %RC, GC
medReadCount=median(temp(:,1)); %median RC
for i=1:101 %looking for median Read counts values to the corresponding GC content
    tempindices=find(temp(:,2)==GCmat(i,1));
    GCmat(i,2)=median(temp(tempindices));
end

for i=1:length(coverage_GC_normalized) % normalization
    positionGC=int32(GCmat(coverage_GC_normalized(i,3)+1,2)); %medReadCount with given GC content
    if positionGC<1 % skip zeros medianReadCount values
        coverage_GC_normalized(i,2)=coverage_GC_normalized(i,2);
    else
        coverage_GC_normalized(i,2)=coverage_GC_normalized(i,2) * ( medReadCount/(positionGC));
    end
end

coverage_GC_normalized=coverage_GC_normalized(:,1:2);
disp('CNV1 GC normalization DONE')

%% BEDGRAPH file, extracting data by 10k rows
bedgraph = tdfread(genmap_file,'tab');
bedgraph = struct2cell(bedgraph);
mappability=[bedgraph{2},bedgraph{3},bedgraph{4}];

%% READ MAPPABILITY NORMALIZATION
% coverage_GC 1sl pozice 2sl RC 3sl GC
coverage_RM=coverage_GC_normalized;
% average RC in mappability count
% mappability 1sl start 2sl stop 3sl MAPscore 4sl average RC
for i=1:length(mappability)
    mappability(i,4)=mean(coverage_RM(mappability(i,1):mappability(i,2),2));
end

% median MAPPscore
medReadCount=median(mappability(:,4)); %median RC for all bins

% MAPPscore matrix
uniqueMAPPscore = unique(mappability(:,3));
MAPPmat=zeros(length(uniqueMAPPscore), 2);
MAPPmat(:,1)=uniqueMAPPscore; %1sl GC%, 2sl median Read Counts with corresponding GC content value
for i=1:length(uniqueMAPPscore)
    value=uniqueMAPPscore(i);
    tempindices=find(mappability(:,3)==value);
    MAPPmat(i,2)=median(mappability(tempindices,4));
end

% mappability 1sl start 2sl stop 3sl MAPscore 4sl average RC 5sl medianRC
% of that mappabilyity score
for i=1:length(mappability)
    findices=find(MAPPmat(:,1)==mappability(i,3));
    mappability(i,5)=MAPPmat(findices,2);
end


%  Normalization
for i=1:length(mappability)
    start=mappability(i,1);
    stop=mappability(i,2);
    coeficient=medReadCount/mappability(i,5);
    coverage_RM(start:stop,2)=coverage_RM(start:stop,2)*coeficient;
    
end
disp('CNV1 Read mappability normalization DONE')

% figure
% plot(coverage_raw(:,2))
% figure
% plot(coverage_GC(:,2))
% figure
% plot(coverage_RM(:,2))
% diff=abs(coverage_GC(:,2)-coverage_RM(:,2));
% figure
% plot(diff)
%% REPLICATION ORIGIN BIAS
coverage_oriC=coverage_RM;

for i=1:length(mappability)
    start=mappability(i,1);
    stop=mappability(i,2);
    distance(i)=abs(mean(start,stop)-oriC_position); %distance from ORIC
    read_count_norm(i)=mean(coverage_oriC(start:stop,2));
    read_count_raw(i)=mean(coverage_raw(start:stop,2));
end

[rho,pval] = corr(distance',read_count_norm','Type','Spearman');
disp(['CNV1 Spearman rank correlation coefficient between distance and read depth: ' num2str(rho)])


% % Fit a poisson distribution
% [b,dev,stats] = glmfit(distance,read_count_norm,'poisson')
% mdl =  fitglm(distance,read_count_norm,'linear','Distribution','poisson')
% mdl =  fitglm(distance,read_count_norm,'Formula','y ~ x1 + x2 + x3')




if rho > 0.9
    %     bp=[1:100000:length(coverage_oriC(:,2))];
    %     coverage_oriC(:,2)=detrend(single(coverage_oriC(:,2)),1,bp,'Continuous',false);
    log2ratio=log2(read_count_raw./read_count_norm);
    mdl =  fitglm(distance,log2(read_count_norm),'y ~ x1 + 1')
    Estimate = mdl.Coefficients.Estimate(2);
    
    for i=1:length(mappability)
        start=mappability(i,1);
        stop=mappability(i,2);
        
        coverage_oriC(start:stop,2)=2.^(log2(single(coverage_oriC(start:stop,2)))+Estimate*single(coverage_oriC(start:stop,2)));
        
    end
else
    
end

disp('CNV1 Replication bias normalization DONE')

%% FINAl output

coverage_final=coverage_oriC(:,1:2);

end