% HISTOGRAM: Displays a histogram of the primary particle sizes.
%
% INPUTS: 
% [STRUCT] = Aggs structure to plot in the histogram
% [INDEX] = Index of aggregate within STRUCT to plot.
%
% OUTPUTS:
% [H] = Histogram struct containing fields used in histogram plot.
%
% Darwin Zhu, 2021-05-26
% ==============================================

function h = agghistogram(struct, index)
dp = struct(index).Pp_manual.dp;

% Sturge's Rule
% numBins = round(1+3.322*log(max(radii))/log(10));

% Rice's Rule
numBins = round(max(dp)^(1/3)*2);

h = histogram(dp, numBins);
title('Primary particle size distribution');
xlabel('Particle dp');
ylabel('Number of particles');
end