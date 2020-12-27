function [filteredSignal,insertSize] = CNV3_pair_reads_distance(bam_file,step)
%% DISTANCE between reads acquired from TLEN column in BAM file
% step - size of window for shifting through bam file (genome coordinates), recommended above 100

InfoStruct = baminfo(bam_file);
ref_ID=InfoStruct.SequenceDictionary.SequenceName;
ref_length=InfoStruct.SequenceDictionary.SequenceLength;

warning('off','all')
% warning
%% DIGGING DATA
filteredSignal=zeros(1,round(ref_length/step));

ii=1;
for k=1:step:round(ref_length/step)*step
    
    BAMStruct=bamread(bam_file,ref_ID,[k k+step]); %BAM BETWEEN COORDINATES
    
    if isempty(BAMStruct) % EMPTY spaces without any reads mapped
        filteredSignal(ii)=0;
    else
        distances=zeros(1,length(BAMStruct),'int32'); %vektor insert size - paired reads distance
        reads_length=zeros(1,length(BAMStruct),'int32'); % vektor delky cteni
        for i=1:length(BAMStruct)
            reads_length(1,i)=length(BAMStruct(i).Sequence);
            distances(1,i)=BAMStruct(i).InsertSize;
        end        
        
        %% Circular genome correction
        % positive values
        distances(distances>ref_length/2)=distances(distances>ref_length/2)-2*reads_length(distances>ref_length/2)-ref_length;
        % negative values
        distances(distances<-ref_length/2)=distances(distances<-ref_length/2)+2*reads_length(distances<-ref_length/2)+ref_length;
        
        %% Outliers
%         idx=isoutlier(single(distances),'median');
%         distances(idx)=[];
        
        %% STATS
        filteredSignal(ii)=mean(abs(distances)); %or median?
        
        %% incrementation of counter
        ii=ii+1;
    end
end

insertSize=mean(abs(filteredSignal));%or median?
filteredSignal=repelem(filteredSignal,step);


disp(['CNV3 Read-pairs Distances computed DONE Insert Size: ' num2str(insertSize)])
%% PLOTTING
figure
plot(filteredSignal,'Color', '#808080')
xlim([-30000 length(filteredSignal)+30000])
ylim([-50 max(filteredSignal)+50])
ylabel('Read-pairs Distance [bp]')
xlabel('Position [Mbp]')
title('Read-pairs Distance')

end

