function [gmean, gstd, ci95] = geomstats(x, w)
% "geomsummary" computes geometric mean, geometric std. dev., and 95% CI
% (log-normal approximation), in weighted or unweighted mode
% -----------------------------------------------------------------------
%
% Inputs:
% x : vector of positive values
% w : (optional) corresponding weights (same length as x, non-negative).
%     If omitted or empty, the estimate is unweighted.
% -----------------------------------------------------------------------
%
% Outputs:
% gmean : (weighted or unweighted) geometric mean
% gstd  : (weighted or unweighted) geometric standard deviation
% ci95  : 95% confidence interval for the geometric mean, [lo hi]
% -----------------------------------------------------------------------

% Ensure input is column vector
x = x(:);

% If no weights provided, compute unweighted statistics
if (nargin < 2) || isempty(w)
    % -------------------- Unweighted case --------------------
    % Sanitize x
    x = x(isfinite(x) & (x > 0));
    n = numel(x);

    if n == 0
        gmean = NaN; gstd = NaN; ci95 = [NaN NaN];
        return
    elseif n == 1
        gmean = x;
        gstd  = 1;         % no dispersion with a single value
        ci95  = [x x];
        return
    end

    y  = log(x);
    mu = mean(y);
    s  = std(y, 0);       % sample std in log-space
    gmean = exp(mu);
    gstd  = exp(s);

    se  = s / sqrt(n);    % standard error of log-mean
    z   = 1.96;
    ci95 = exp(mu + z * [-1 1] * se);

else
    % -------------------- Weighted case ----------------------
    w = w(:);

    % Joint mask on x and w
    ok = isfinite(x) & (x > 0) & isfinite(w) & (w > 0);
    x = x(ok);
    w = w(ok);

    if isempty(x)
        gmean = NaN; gstd = NaN; ci95 = [NaN NaN];
        return
    end

    % Normalize weights
    w = w / sum(w);

    % Weighted geometric mean
    y = log(x);
    mu = sum(w .* y);
    gmean = exp(mu);

    % Weighted geometric standard deviation
    variance_log = sum(w .* (y - mu).^2);
    gstd = exp(sqrt(variance_log));

    % Kish effective sample size
    n_eff = 1 / sum(w.^2);

    if n_eff > 1
        se  = sqrt(variance_log) / sqrt(n_eff);
        z   = 1.96;
        ci95 = exp(mu + z * [-1 1] * se);
    else
        ci95 = [gmean gmean];
        warning('gmci:LowEffectiveN', ...
            'Effective sample size n_eff = %.2f ≤ 1. CI is not reliable.', n_eff);
    end
end
end
