
% VIZ_DADP  Generate a formatted plot of da versus dp. 
% This includes a fit relation and the universal relation of Olfert and Rogak.
% Author: Timothy Sipkens
%=========================================================================%

function [] = viz_dadp(Aggs_da, dp)

% Parse inputs
if isstruct(Aggs_da); da = [Aggs_da.da]; dp = [Aggs_da.dp];
else da = Aggs_da; end


% Plot data
loglog(da, dp, '.', 'Color', [0.12,0.59,0.96]);
hold on;


% Get current figure limits
xlims = xlim; ylims = ylim;
t0 = min([xlims(1),ylims(1)]);
t1 = max([xlims(2),ylims(2)]);


% Plot dp-da relation
p1 = polyfit(log10(da), log10(dp),1);
p1_val = 10 .^ polyval(p1, log10([t0,t1]));
loglog([t0,t1], p1_val, 'Color', [0.12,0.59,0.96]);


% Plot universal relation from Olfert and Rogak
loglog([t0,t1], 10.^(log10(17.8) + 0.35.*log10([t0,t1]./100)), 'k--');


%-- Plot 2-sigma ellipse -------------------------------------------------%
mu = [mean(log10(da)), mean(log10(dp))]; % mean for ellipse center
Sigma = cov(log10(da), log10(dp)); % covariance of plotting ellipse

[V, D] = eig(Sigma.*2);
t = linspace(0, 2*pi); % points around the edge of the circle
a = (V*sqrt(D)) * [cos(t(:))'; sin(t(:))']; % points on the ellipse

loglog(10.^(a(1, :)+mu(1)), 10.^(a(2, :)+mu(2)), '--', ...
    'Color', [0.12,0.59,0.96]);
%-------------------------------------------------------------------------%


% Plot polygon of off-limit particles (i.e., where da > dp)
% Edge of this region corresponds to single primary particle aggregates.
xlims = xlim; ylims = ylim;
fill([xlims(1), xlims(1), ylims(2)], ...
    [xlims(1), ylims(2), ylims(2)], ...
    [0.95,0.95,0.95], 'EdgeColor', [0.5,0.5,0.5]);


% Rearrange order, so shaded region is in background
h = get(gca, 'Children');
set(gca, 'Children', [h(2:end); h(1)]);

xlabel('d_a [nm]');
ylabel('d_p [nm]');
legend({'d_p > d_a', 'Data', 'Power law fit', ...
    'Universal relation', '2-sigma ellipse'}, ...
    'Location', 'southeast');


hold off;
if nargout==0; clear h; end % remove output if none requested

end

