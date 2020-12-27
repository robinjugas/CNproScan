function [outliers, idx] = gesd(coverageSignal, estimatedOutliers)
% GESD Two-sided Generalized (extreme Studentized deviate) ESD test
% coverageSignal vector of coverage
% estimatedOutliers number of estimated over threshold values

%% COMPUTATION
alpha = 0.05;
R=zeros(1,estimatedOutliers);
Lambda=zeros(1,estimatedOutliers);
outliers=zeros(1,estimatedOutliers);

coverageSignal=double(coverageSignal);

n = length(coverageSignal);
xCopy = coverageSignal;
outliersNumber=0;
for i = 1 : estimatedOutliers
    xMean = mean(coverageSignal);
    xSD = std(coverageSignal);
    
    R(i) = max( abs(coverageSignal - xMean) / xSD );
    
    % compute critical value
    p = 1 - alpha /2 /(n - i + 1);
    t = tinv(p, n - i - 1);
    Lambda(i) = t * (n - i) / sqrt((n - i - 1 + t^2) * (n - i + 1));
    
    if R(i) > Lambda(i)
        outliersNumber = i;
    end
    
    [~, idxMax] = max( abs(coverageSignal - xMean) );
    outliers(i) = coverageSignal(idxMax);
    coverageSignal(idxMax) = [];
end

idx = zeros(1, outliersNumber);

for i = 1 : outliersNumber
    idx(i) = find(xCopy == outliers(i), 1);
    xCopy( idx(i) ) = NaN;
end

if outliersNumber > 0
    outliers = outliers(1 : outliersNumber);
    [idx,order]=sort(idx);
    outliers=outliers(order);
else
    outliers=[];
    idx=[];
end




end