function [gmean, gstd] = weighted_geomean(x, w)
% weighted_geomean computes weighted geometric mean and geometric std. dev.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   x : vector of positive values
%   w : corresponding weights (same length as x, non-negative)
% ----------------------------------------------------------------------- %
%
% Outputs:
%   gmean : weighted geometric mean
%   gstd  : weighted geometric standard deviation
% ----------------------------------------------------------------------- %

% Ensure inputs are column vectors
x = x(:);
w = w(:);

% Normalize weights to sum to 1
w = w / sum(w);

% Weighted geometric mean
gmean = exp(sum(w .* log(x)));

% Weighted geometric standard deviation
variance_log = sum(w .* (log(x) - log(gmean)).^2);
gstd = exp(sqrt(variance_log));

end