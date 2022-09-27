function pp = manu_pp_0(fnam0, da, fad0, opts)
% "manu_pp" calculates the properties of primary particles extracted from...
%   ...manual sizing of a real TEM aggregate.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   fnam0: A cell array of name of the files to be imported
%   da: Projected area diameter of aggs
%   fad0: folder address (e.g., 'inputs\manual\')
%   n_hyb: Number of hybridity regions
%   opts: Data import options
% ----------------------------------------------------------------------- %
%
% Output:
%   pp: Primary particle data structure
% ----------------------------------------------------------------------- %

n_agg = length(fnam0); % Number of aggs to be analyzed

% Initialize the pp structure
pp = struct();
pp.mu_d = zeros(n_agg,2); % Mean primary particle diameter of aggs (regular + geometric)
pp.std_d = zeros(n_agg,2); % Standard deviation of pp size (regular + geometric)
pp.n = zeros(n_agg,1); % Number of pps
pp.d = cell(n_agg,1); % Projected area diameter of pps
pp.r = cell(n_agg,1); % pp centroid locations
pp.c = cell(n_agg,1); % pp circularities

% Set defaults for dp vs da plotting
if ~exist('opts', 'var') 
    opts = struct();
end

if ~isfield(opts, 'visual')
    opts.visual = [];
end
opts_visual = opts.visual;
if isempty(opts_visual)
    opts_visual = 'off'; % Default not to plot the outputs
end

if ~isfield(opts, 'lbl')
    opts.lbl = cell(1, n_agg);
end

% if ~isfield(opts, 'nhyb')
%     opts.nhyb = []; % number of hybridity regions
% end

% % set colormap
% if ~isempty(nhyb)
%     mc = colormap(turbo);
%     ii = round(1 + (length(mc) - 1) .*...
%         (0.05 : 0.9 / (length(unique(n_hyb)) - 1) : 0.95)');
%     mc = mc(ii,:);
%     mc = flip(mc,1);
%     
%     [nhyb_s, I_s] = sort(nhyb);
%     nhyb_us = unique(nhyb);
%     
% end
% 
% % adjust coloring order based on number of hybridity regions

pp.fname = fnam0; % store file names

% Import the data and calculate properties
for i = 1 : n_agg
    fad = strcat(fad0, fnam0{i}, '.csv'); % Making file address
    ppdat = readmatrix(fad); % Scanning data
    
    pp.n(i) = size(ppdat, 1) - 2; % Removing the header and the scale detector 
    pp.d{i} = 2 * sqrt(ppdat(3:end, 2) / pi);
    pp.r{i} = ppdat(3:end, 3:4);
    pp.c{i} = 2 * sqrt(pi * ppdat(3:end, 2)) ./ ppdat(3:end, 5);
    
    pp.mu_d(i,:) = [mean(pp.d{i}), geomean(pp.d{i})];
    pp.std_d(i,:) = [std(pp.d{i}), GEOSTD(pp.d{i})];
end

if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 600, 600])
        h.Position = [0, 0, 600, 600];
    end
    set(h, 'color', 'white')
    
    % Plot universal correlation
    dpp_uc = linspace(5, 65, 1000);
    da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);
    p1 = plot(da_uc, dpp_uc, 'Color', [0.5 0.5 0.5],...
        'LineStyle', '-.', 'LineWidth', 2.5);
    hold on
    
    % Plot manually sized aggs
%     if ~isempty(nhyb)
%         iii = find(nhyb == nhyb_us(i));
%     else
        p2 = scatter(da, pp.mu_d(:,2), 25, [0.8500 0.3250 0.0980], '^');
    
    if (~isfield(opts, 'eb')) || isempty(opts.eb)
        opts.eb = 'on'; % default to plot errorbars
    end
    opts_eb = opts.eb;
    if ismember(opts_eb, {'ON', 'On', 'on'})
        e_p = pp.mu_d(:,2) .* abs(pp.std_d(:,2) - 1);
        e_n = pp.mu_d(:,2) .* abs(1 - 1 ./ pp.std_d(:,2));
        eb = errorbar(da, pp.mu_d(:,2), e_n, e_p, 'Marker', 'none',...
            'LineStyle', 'none');
        eb.Color = [0.8500 0.3250 0.0980];
    end
    
    if (~isfield(opts, 'label')) || isempty(opts.label)
        opts.label = 'on'; % default to plot std labels
    end
    opts_label = opts.label;
    if ismember(opts_label, {'ON', 'On', 'on'})
        text(da, pp.mu_d(:,2), num2str(pp.std_d(:,2), '%.2f'),...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
            'FontSize', 8)
        text(da, pp.mu_d(:,2), opts.lbl, 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'right', 'FontSize', 8)        
    end        

    box on
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    xlabel('d_a (nm)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
        'FontSize', 14)
    ylabel('$\overline {d_{pp}} (nm)$', 'Interpreter', 'latex',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 14)
    xlim([5, inf])
    ylim([5, 65])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    legend([p1, p2], {'Universal correlation', 'Manual'},...
        'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
    title('Primary particle size vs projected area equivalent size',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
end

end


function gsd = GEOSTD(x, flag, dim)
% For vectors, GEOSTD(X) calculates the geometric standard deviation of
% input X.  For  matrices GEOSTD calculates the geometric standard
% deviation over each column of X.  For N-D arrays, GEOSTD calculates the
% geometric standard deviation over the first non-singleton dimension of X.
% 
% 
% WARNING: By default GEOSTD normalises by N, where N is the sample size.
% This is different to the way the built-in MATLAB function STD operates.
% 
% 
% GEOSTD(X, 0) calculates the geometric standard deviation as above, with
% normalisation factor (N-1).  GEOSTD(X, 1) works like GEOSTD(X).
% 
% 
% GEOSTD(X, [], DIM) calculates the geometric standard deviation along
% dimension DIM of X.  Use a flag of 0 to normalise by factor (N-1), or a
% flag of 1 to normalise by factor N.
% 
% 
% NOTE: Class type checking and error handling are conducted within EXP,
% STD and LOG.
%
% 
% EXAMPLE:  X = 10*rand(5);
%                      geostd(X)
%                      ans =
%                     
%                         1.1858    1.8815    1.8029    4.1804    2.5704
% 
% 
%   Class support for input X:
%      float: double, single
% 
% 
%   See also GEOMEAN (stats toolbox), STD.
% 
% 
% $ Author: Richie Cotton $     $ Date: 2006/03/17 $


% Basic error checking of inputs
if nargin < 1
    error('geostd:notEnoughInputs', 'This function requires at least one input.');
elseif any(x(:) < 0)
    error(geostd:badData', 'All data values must be positive.');
end

% Setup default flag where required
if nargin < 2 || isempty(flag)
    flag = 1;
end

% If dimension is not specified, find first non-singleton dimension
if nargin < 3 || isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

% Turn off warnings regarding log of zero, since this is an artifact of the
% technique use to calculate the gsd
lozwarning = warning('off', 'MATLAB:log:logOfZero');

% Calculate geometric std dev using 
% "log of geometric std dev of data = arithmetic std dev of log of data"
gsd = exp(std(log(x), flag, dim));

% Reset warning value
warning(lozwarning);
end