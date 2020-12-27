function [indicationHigher,indicationLower] = CNV4_pair_reads_distance_thresholding(signalDistance,insertSize)
%READS_DISTANCE_PROCESSING extension of areas AND Thresholding to create
%binary vector representing areas with higher reads distance
%threshold = integer or 'statistics' with raw_data computed or empty
insertSize=abs(insertSize);
signalDistance=abs(signalDistance);
upper_threshold=insertSize+2*std(signalDistance);
lower_threshold=insertSize-2*std(signalDistance);
extension=100;

%% THRESHOLDING HIGHER
[~,idxHigh]=find(signalDistance >= upper_threshold);
valHigh=signalDistance(idxHigh);
%binary
indicationHigher=zeros(1,length(signalDistance));
indicationHigher(signalDistance >= upper_threshold) = 1;
%extension
start=strfind([0 indicationHigher],[0 1]); %najde prechod 0 1
stop=strfind([indicationHigher 0],[1 0]); %najde prechod 1 0
for i=1:length(start)
    indicationHigher(start(i)-extension:start(i))=ones(1,extension+1);
    indicationHigher(stop(i):stop(i)+extension)=ones(1,extension+1);
end
% figure;plot(indicationHigher)

%% THRESHOLDING LOWER
[~,idxLow]=find(signalDistance < lower_threshold);
valLow=signalDistance(idxLow);
%binary
indicationLower=zeros(1,length(signalDistance));
indicationLower(signalDistance < lower_threshold) = 1;
%extension
start=strfind([0 indicationLower],[0 1]); %najde prechod 0 1
stop=strfind([indicationLower 0],[1 0]); %najde prechod 1 0
for i=1:length(start)
    indicationLower(start(i)-extension:start(i))=ones(1,extension+1);
    indicationLower(stop(i):stop(i)+extension)=ones(1,extension+1);
end
% figure;plot(indicationLower)
%% PLOTTING POLISHED PEAKS
figure
plot(signalDistance,'Color', '#808080')
hold on
plot(idxHigh,valHigh,'xb')
hold on
plot(idxLow,valLow,'xr')

xlim([-30000 length(signalDistance)+30000])
ylim([-50 max(signalDistance)+50])

ylabel('Read-pairs Distance [bp]')
xlabel('Position [Mbp]')
title('Read-pairs Distance')
% title(['Pair-read Distance - Insert Size: ' num2str(insertSize)])
legend('Read-pairs Distance','High Outliers','Low Outliers')
legend('boxoff')
box off

%%
disp('CNV4 Read-pairs Distance thresholding DONE')

end

