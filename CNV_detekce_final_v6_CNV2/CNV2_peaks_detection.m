function [peaksPolished,indicationSignal] = CNV2_peaks_detection(coverageSignal)
%%Peak detection
% coverageSignal vector of coverage data
disp('Version 6.0')

% SETTINGS
step=11; % STEP size, how many values are taken into slope equation
peakDistanceThreshold=20; %INTEGER, distance between peaks threshold, if
% shorter, peaks are merged
cutoff=500; % length of segments in the start and end to be replaced - fake deletions

%% DELETING crossovers at start and end

coverageSignal(1:cutoff)=repelem(mean(coverageSignal),length(coverageSignal(1:cutoff)));
coverageSignal(end-cutoff:end)=repelem(mean(coverageSignal),length(coverageSignal(end-cutoff:end)));

%% ZERO COVERAGE omitted and detected
modifiedSignal=double(coverageSignal);
idxDeletions=find(modifiedSignal==0);
modifiedSignal(idxDeletions)=mean(modifiedSignal);


%% ESTIMATING NUMBER  OF OUTLIERS - MODIFIED Z-SCORE
[outliersZSCORE, idxZSCORE] = mzscore(modifiedSignal);
estimatedNumOutliers=length(idxZSCORE);

%% thresholding of MZSCORE by three-sigma rule to reduce number of possible
% outliers
% if estimatedNumOutliers>20000
%     sigmaruleP=mean(modifiedSignal)+3*std(modifiedSignal);
%     sigmaruleN=mean(modifiedSignal)-3*std(modifiedSignal);
%     % idxZSCORE_orig=idxZSCORE;outliersZSCORE_orig=outliersZSCORE;
%     idxZSCORE(outliersZSCORE<=sigmaruleP & outliersZSCORE>=sigmaruleN)=[];
%     outliersZSCORE(outliersZSCORE<=sigmaruleP & outliersZSCORE>=sigmaruleN)=[];
%     estimatedNumOutliers=length(idxZSCORE);
% end

disp(['CNV2 M-Zscore estimated outliers: ' num2str(estimatedNumOutliers)])

%% GESD Generalized Extreme Studentized Deviation (GESD)
[outliersGESD, idxGESD]=gesd(modifiedSignal,estimatedNumOutliers);

if isempty(outliersGESD) || isempty(idxGESD)
    peaksPolished=[];
    indicationSignal=0;
    disp('CNV2 No outliers detected')
    figure
    plot(coverageSignal)    
    ylabel('Coverage')
    xlabel('Position [bp]')
    title('Coverage - no CNV detected') 
    box off
    return
end

disp('CNV2 GESD outliers detection DONE')

%% ADDING POTENTIALLY SKIPPED DELETIONS
% idx=idxGESD;
% deletions=find(coverageSignal==0);
idx = union(idxGESD,idxDeletions);

%% MERGING CLOSE PEAKS TOGETHER

differences=diff(idx);
kk=find(differences>1 & differences<peakDistanceThreshold);
for i=1:length(kk)
    pom=idx(kk(i)+1)-idx(kk(i));
    vec=idx(kk(i))+1 : idx(kk(i))+pom-1;
    
    idx=[ idx(1:kk(i)) vec  idx(kk(i)+1:end) ];
    
    korekce=pom-1; %korekce - sekvence o korekci delsi, musim odecist
    kk=kk+korekce;
end

%% PREPROCESSING of peaks into cell datatype
outliers=coverageSignal(idx);

peaks=cell(0,0);
differences=diff(idx);
k=find(differences>1);
for i=1:length(k)+1
    if i==1 %FIRST peak
        peaks{i}(:,1)=idx(1 : k(i));
        peaks{i}(:,2)=outliers(1 : k(i));
    elseif i==length(k)+1 %LAST peak
        peaks{i}(:,1)=idx(k(i-1)+1 : end);
        peaks{i}(:,2)=outliers(k(i-1)+1 : end);
    else %IN BETWEEN peaks
        peaks{i}(:,1)=idx( k(i-1)+1: k(i));
        peaks{i}(:,2)=outliers( k(i-1)+1: k(i));
    end
end


%% PEAK BORDERS DETECTION based on slope of a peak - works only on POSITIVE PEAKS, NOT NEGATIVE PEAKS - DELETIONS
peaksPolished=cell(0,0);
avg_signal=mean(coverageSignal);%average coverage of genome

for i=1:length(peaks)
    start=peaks{i}(1,1);
    stop=peaks{i}(end,1);
    avg_cov_peak=mean(coverageSignal(start:stop));%CNVcoverage
    
    if avg_cov_peak<=avg_signal %DELETIONS \/
        %SLOPE \ peak start
        count1=0;
        trigger1=0;
        while trigger1<5
            if start<=1 || start-count1*step-step<=1
                break
            else
                vektor=coverageSignal(start-count1*step-step:start-count1*step);
                mod_vektor=mode(vektor);
            end
            if find(~vektor) %if zero in vector, ends - cause zero coverage
                break
            end
            
            count1 = count1 + 1;
            if mod_vektor>=avg_cov_peak
                trigger1=trigger1+1;
            end
        end
        % SLOPE / peak stop
        count2=0;
        trigger2=0;
        while trigger2<5
            if stop>=length(coverageSignal) || stop+count2*step+step>=length(coverageSignal)
                break
            else
                vektor=coverageSignal(stop+count2*step : stop+count2*step+step);
                mod_vektor=mode(vektor);
            end
            
            if find(~vektor) %if zero in vector, ends - cause zero coverage
                break
            end
            
            count2 = count2 + 1;
            if mod_vektor>=avg_cov_peak
                trigger2=trigger2+1;
            end
        end
        vektor1=start-count1*step : start-1;
        vektor=peaks{i}(:,1)';
        vektor2=stop+1:stop+count2*step;
        peaksPolished{i}(:,1)=[vektor1 vektor vektor2]';
        peaksPolished{i}(:,2)=coverageSignal(peaksPolished{i}(:,1))';
        
    else %NOT DELETIONS /\
        %SLOPE UPWARDS peak start /
        count1=0;
        trigger1=0;
        while trigger1<5
            if start<=1 || start-count1*step-step<=1
                break
            else
                vektor=coverageSignal(start-count1*step-step:start-count1*step);
            end
            
            if find(~vektor) %if zero in vector, ends - cause zero coverage
                break
            end
            slope=(vektor(1)-vektor(end))/-length(vektor);
            count1 = count1 + 1;
            if slope<=0
                trigger1=trigger1+1;
            end
        end
        % SLOPE DOWNWARDS peak stop \
        count2=0;
        trigger2=0;
        while trigger2<5
            if stop>=length(coverageSignal) || stop+count2*step+step>=length(coverageSignal)
                break
            else
                vektor=coverageSignal(stop+count2*step : stop+count2*step+step);
            end
            
            if find(~vektor) %if zero in vector, ends - cause zero coverage
                break
            end
            slope=(vektor(1)-vektor(end))/-length(vektor);
            count2 = count2 + 1;
            if slope>=0
                trigger2=trigger2+1;
            end
        end
        vektor1=start-count1*step : start-1;
        vektor=peaks{i}(:,1)';
        vektor2=stop+1:stop+count2*step;
        peaksPolished{i}(:,1)=[vektor1 vektor vektor2]';
        peaksPolished{i}(:,2)=coverageSignal(peaksPolished{i}(:,1))';
    end
    
end

%% BINARY indication vector, 1 where is the peak present
indicationSignal=zeros(1,length(coverageSignal));
for i=1:length(peaksPolished)
    indicationSignal(peaksPolished{i}(:,1)')=1;
end

disp('CNV2 Peaks detection DONE')

%% PLOTTING Zscore vs GESD
figure
plot(coverageSignal,'Color', '#808080')
hold on
plot(idxZSCORE,coverageSignal(idxZSCORE),'xb')
hold on
plot(idxGESD,coverageSignal(idxGESD),'xg')
hold on
plot(idxDeletions,coverageSignal(idxDeletions),'xr')

ylabel('Coverage')
xlabel('Position [Mbp]')
xlim([-30000 length(coverageSignal)+30000])
ylim([0 max(coverageSignal)+50])
title('Outliers in Coverage Signal')
box off
legend('Coverage','M-Zscore Outliers','GESD Outliers','Zero Coverage')
legend('boxoff')

%% PLOTTING POLISHED PEAKS
figure
plot(coverageSignal,'Color', '#808080')
hold on
for i=1:length(peaksPolished)    
    plot(peaksPolished{i}(:,1),peaksPolished{i}(:,2),'b')
    hold on
end
xlim([-30000 length(coverageSignal)+30000])
ylim([0 max(coverageSignal)+50])
ylabel('Coverage')
xlabel('Position [Mbp]')
title('CNV Regions')
box off
legend('Coverage','CNV Peaks')
legend('boxoff')

% %% NORM DIST PLOTS
% withoutOutliersGESD=coverageSignal;
% withoutOutliersGESD(idxGESD)=[];
% withoutOutliersDeletions=coverageSignal;
% withoutOutliersDeletions(idxDeletions)=[];
% withoutOutliersCNV=coverageSignal;
% withoutOutliersCNV(indicationSignal==1)=[];
% 
% figure
% subplot(2,2,1)
% normplot(coverageSignal)
% title('Original Coverage Distribution')
% xlabel('Coverage')
% 
% subplot(2,2,2)
% normplot(withoutOutliersCNV)
% title('Coverage Distribution without CNV regions')
% xlabel('Coverage')
% 
% subplot(2,2,3)
% normplot(withoutOutliersDeletions)
% title('Coverage distribution without zero values' )
% xlabel('Coverage')
% 
% subplot(2,2,4)
% normplot(withoutOutliersGESD)
% title('Coverage distribution without GESD outliers')
% xlabel('Coverage')


end

