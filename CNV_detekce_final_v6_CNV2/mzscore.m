function [outliers, idx] = mzscore(x)
% MZSCORE Modified Z-score to screen data for outliers by Iglewicz and Hoaglin
% Mi=0.6745(xi?x~)MAD

%% check of input arguments

if nargin == 0
    error('Input arguments are NOT passed to function')
end

% number of input arguments must be 1
narginchk(1, 1)

%%  Iglewicz and Hoaglin outlier test  Z-score Calculation

n = length(x);
M = zeros(1, n);

medianX = median(x);
MAD = median( abs(x - medianX) );

for i = 1 : n
    M(i) = 0.6745 * abs(x(i) - medianX) / MAD;
end

idx = find(M > 3.5);
outliers = x(idx);

end