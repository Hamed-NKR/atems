% HISTOGRAM: Displays a histogram of the primary particle sizes.
% Darwin Zhu, 2021-05-26
% ==============================================

function h = histogram(struct, index)
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