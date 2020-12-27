function [CNVtable,CNVseq] = CNV5_detection_output(peaksPolished,binaryPeaks,binaryDistanceHigher,binaryDistanceLower,coverageSignal,FASTAfile)
%% CNV_DETECTION Summary of this function goes here
[~,ref_seq]=fastaread(FASTAfile); %length of reference used to BWA

%%
CNVseq=struct('Sequence',{},'Header',{});
%% DELETING crossovers at start and end
cutoff=500;
coverageSignal(1:cutoff)=repelem(mean(coverageSignal),length(coverageSignal(1:cutoff)));
coverageSignal(end-cutoff:end)=repelem(mean(coverageSignal),length(coverageSignal(end-cutoff:end)));

avg_signal=mean(coverageSignal);%average coverage of genome


%% TABLE GENERATING
CNVtable = cell(1,10);
CNVtable{1}='ID';CNVtable{2}='StartPosition';CNVtable{3}='StopPosition';CNVtable{4}='CNVlength';
CNVtable{5}='CNVcoverage';CNVtable{6}='CNVrelativeCoverage';CNVtable{7}='CNVtype';CNVtable{8}='Sequence';
CNVtable{9}='GenomeAverageCoverage:';CNVtable{10}=avg_signal;

%% LENGTH DIFFERENCE TESTING
dif=length(binaryPeaks)-length(binaryDistanceHigher);
if dif>0
    binaryDistanceHigher=[binaryDistanceHigher zeros(1,dif)];
    binaryDistanceLower=[binaryDistanceLower zeros(1,dif)];
elseif dif<0
    binaryPeaks=[binaryPeaks zeros(1,-dif)];
end
%% IDX of 1 in vector
[~,distanceHigh_idx]=find(binaryDistanceHigher);% elevated read distance coordinates
[~,distanceLow_idx]=find(binaryDistanceLower);% elevated read distance coordinates

%% SORTING & PLOTTING
figure
plot(coverageSignal,'Color', '#808080')
hold on
for i=1:length(peaksPolished)
    x=peaksPolished{1,i}(:,1);
    y=peaksPolished{1,i}(:,2);
    avg_cov_peak=mean(y);%CNVcoverage
    rel_cov_peak=(avg_cov_peak/avg_signal)*100;%CNVrelativeCoverage
    start=x(1);
    stop=x(end);
    peak_length=length(x);
    varHigh=ismember(x,distanceHigh_idx);%in area of distorted reads Higher Insert Size
    varLow=ismember(x,distanceLow_idx);%in area of distorted reads lower Insert Size
    % type of CNV
    if avg_cov_peak<=avg_signal && ~isempty(find(varHigh,1)) %LOW COVERAGE+HIGH DISTANCE SAME IN PLOT
        type='Deletion - zero coverage+discordant reads';
        hold on
        plot(x,y,'c')
    elseif avg_cov_peak<=avg_signal %LOW COVERAGE only
        type='Deletion - zero coverage';
        hold on
        plot(x,y,'xr')
    elseif avg_cov_peak>=avg_signal && ~isempty(find(varHigh,1)) %HIGH COVERAGE+HIGH DISTANCE
        type='Interspersed CNV';
        hold on
        plot(x,y,'g')
    else
        type='Tandem CNV';
        hold on
        plot(x,y,'b')
    end
    
    % table
    seq(i).Sequence=ref_seq( start : stop );
    ID=[num2str(i) char(randi([65,90])) char(randi(+'AZ'))];
    seq(i).Header=[ 'CNVid:' ID ' ReferenceSequencePostion: ' num2str(start) ' : ' num2str(stop) ];
    % +1 because of header row
    CNVtable{i+1,1}= ID; %ID
    CNVtable{i+1,2}=start; %START
    CNVtable{i+1,3}=stop; %STOP
    CNVtable{i+1,4}=peak_length;%LENGTH
    CNVtable{i+1,5}=avg_cov_peak;%CNVcoverage
    CNVtable{i+1,6}=rel_cov_peak;%CNVrelativeCoverage
    CNVtable{i+1,7}=type;
    CNVtable{i+1,8}=seq(i).Sequence;
    
    
end
%% PLOTTING
xlim([-30000 length(coverageSignal)+30000])
ylim([0 max(coverageSignal)+50])

ylabel('Coverage')
xlabel('Position [Mbp]')
title('Detected CNVs')
box off


annotation('textbox',[.13 .83 .1 .1],'String','Deletion+DR','FitBoxToText','on','Color','y','EdgeColor','none');
annotation('textbox',[.13 .80 .1 .1],'String','Deletion','FitBoxToText','on','Color','r','EdgeColor','none');
annotation('textbox',[.13 .77 .1 .1],'String','Interspersed','FitBoxToText','on','Color','g','EdgeColor','none');
annotation('textbox',[.13 .74 .1 .1],'String','Tandem','FitBoxToText','on','Color','b','EdgeColor','none');

disp('CNV5 Results DONE')

end
