function [Aggs_freq, Aggs_bin, h] = freqdist(Aggs, n_bin)
% FREQDIST plots frequncy distribution plots of aggregates based on their
% different morphological properties (i.e., optical depth, circularity,
% and size). 
% ----------------------------------------------------------------------- %
% Aggs: Aggregates data structure
% n_bin: number of post-processing bins (3 member vector for size,...
%   ...circularity and optical depth)
% Aggs_freq: frequency of 1d and 2d size, circularity and optical depth
%   ...bins
% Aggs_bin: bin labels of each aggregate
% h: output figure handle
% ----------------------------------------------------------------------- %

n_agg = length(Aggs.id); % total number of aggregates

% compile the properties from the aggregate population
da = cat(1, Aggs.da); % area equivalent diameter
ci = cat(1, Aggs.circularity); % circularity
od = cat(1, Aggs.zbar_opt); % optical depth

% range of properties
del_da = [min(da), max(da)];
del_ci = [min(ci), max(ci)];
del_od = [min(od), max(od)];

% Assigning the number of bins if not given
if ~exist('n_bin', 'var') || isempty(n_bin)
    n_bin = [5, 5, 5];
elseif length(n_bin) == 1
    n_bin = repmat(n_bin, 1, 3);
elseif length(n_bin(:)) ~= 3
    error('Invalid number of bins! (should be a 3 member vector')
end

% discretize the domains of post-processing variables
da_bin = del_da(1) + (del_da(2) - del_da(1)) .* (0 : 1 / n_bin(1) : 1);
ci_bin = del_ci(1) + (del_ci(2) - del_ci(1)) .* (0 : 1 / n_bin(2) : 1);
od_bin = del_od(1) + (del_od(2) - del_od(1)) .* (0 : 1 / n_bin(3) : 1);

% initialize the bin centers
da_c = zeros(n_bin(1),1);
ci_c = zeros(n_bin(2),1);
od_c = zeros(n_bin(3),1);

% initialize the frequency of occurrence arrays
Aggs_freq = {zeros(n_bin(1),1), zeros(n_bin(2),1), zeros(n_bin(3),1);...
    zeros(n_bin(1), n_bin(2)), zeros(n_bin(1), n_bin(3)),...
    zeros(n_bin(2), n_bin(3))};

% initialize the array to label the aggregates for different bin sets 
Aggs_bin = zeros(n_agg,3);


%%% assign the aggregates to the bins & count for 1d frquencies

% da_1d
for i = 1 : n_bin(1)
    if i == 1 % The first bin (only upper limit)
        ii = da < da_bin(i+1);
        da_c(i) = (da_bin(i) + da_bin(i+1)) / 2;
    elseif i < n_bin(1)
        ii = (da >= da_bin(i)) && (da < da_bin(i+1));
        da_c(i) = (da_bin(i) + da_bin(i+1)) / 2;
    else % the last bin (all the remaining)
        ii = Aggs_bin(:,1) == 0;
    end
    
    Aggs_freq{1,1}(i) = 100 * nnz(ii) / n_agg; % binned size freq. (%)
    Aggs_bin(ii,1) = i; % column 1 to be size bin labels
end

% ci_1d
for j = 1 : n_bin(2)
    if j == 1
        jj = ci < ci_bin(j+1);
        ci_c(j) = (ci_bin(j) + ci_bin(j+1)) / 2;
    elseif j < n_bin(2)
        jj = (ci >= ci_bin(j)) && (ci < ci_bin(j+1));
        ci_c(j) = (ci_bin(j) + ci_bin(j+1)) / 2;
    else
        jj = Aggs_bin(:,2) == 0;
    end
    
    Aggs_freq{1,2}(j) = 100 * nnz(jj) / n_agg; % circularity freq. (%)
    Aggs_bin(jj,2) = j; % column 2 to be circularity bin labels    
end

% od_1d
for k = 1 : n_bin(3)
    if k == 1
        kk = od < od_bin(k+1);
        od_c(k) = (od_bin(k) + od_bin(k+1)) / 2;
    elseif k < n_bin(3)
        kk = (od >= od_bin(k)) && (od < od_bin(k+1));
        od_c(k) = (od_bin(k) + od_bin(k+1)) / 2;
    else
        kk = Aggs_bin(:,3) == 0;
    end
    
    Aggs_freq{1,3}(k) = 100 * nnz(kk) / n_agg; % optical depth freq. (%)
    Aggs_bin(kk,3) = k; % column 3 to be optical depth bin labels
end

%%%


%%% identify 2d frequencies

% da & ci
for i = 1 : n_bin(1)
    for j = 1 : n_bin(2)
        Aggs_freq{2,1}(i,j) = 100 * nnz((Aggs_bin(:,1) == i) &...
            (Aggs_bin(:,2) == j)) / n_agg;
    end
end

% da & od
for i = 1 : n_bin(1)
    for k = 1 : n_bin(3)
        Aggs_freq{2,2}(i,k) = 100 * nnz((Aggs_bin(:,1) == i) &...
            (Aggs_bin(:,3) == k)) / n_agg;
    end
end

% ci & od
for j = 1 : n_bin(2)
    for k = 1 : n_bin(3)
        Aggs_freq{2,3}(j,k) = 100 * nnz((Aggs_bin(:,2) == j) &...
            (Aggs_bin(:,3) == k)) / n_agg;
    end
end

%%%


%%% plotting results

figure;
h = gcf;
h.Position = [0, 0, 900 * 2, 900]; % Position and size
set(h, 'color', 'white'); % Background color

% set the figure layout 
tt = tiledlayout(3, 3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% Colormaps for different variables
cm1 = winter; % for da
ll1 = round(10 + (length(cm1) - 10) .* (0 : 1 / n_bin(1) : 1)');
cm1 = cm1(ll1,:);

cm2 = autumn; % for ci
ll2 = round(10 + (length(cm2) - 10) .* (0 : 1 / n_bin(2) : 1)');
cm2 = cm2(ll2,:);

cm3 = gray; % for od
ll3 = round(10 + (length(cm3) - 10) .* (0 : 1 / n_bin(3) : 1)');
cm3 = cm3(ll3,:);

% a. the 1d frequency distribution barcharts

% a.1. nbar vs da
nexttile
bc1 = bar(da_c, Aggs_freq{1,1}, 'FaceColor', 'flat');
bc1.CData = cm1;
axis equal
box on
xlim(del_da)
ylim([0 1])
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel({'$\bar{n}$', ' (%)'}, 'Interpreter', 'Latex',...
    'FontName', 'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')

% a.2. nbar vs ci
nexttile
bc2 = bar(ci_c, Aggs_freq{1,2}, 'FaceColor', 'flat');
bc2.CData = cm2;
axis equal
box on
xlim(del_ci)
ylim([0 1])
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('$\bar{n}$', 'Interpreter', 'Latex', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')

% a.3. nbar vs od
nexttile
bc3 = bar(od_c, Aggs_freq{1,3}, 'FaceColor', 'flat');
bc3.CData = cm3;
axis equal
box on
xlim(del_od)
ylim([0 1])
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('z_o_p_t (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('$\bar{n}$', 'Interpreter', 'Latex', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')

% b. the 2d frequency distribution heatmaps

% b.1. nbar vs da & ci
nexttile
hm1 = imagesc(da_c, ci_c, (Aggs_freq{2,1})');
axis equal
box on
xlim(del_da)
ylim(del_ci)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb1 = colorbar;
cb1.Label.String = '$\bar{n}$';
cb1.Label.Interpreter = 'Latex';
cb1.Label.FontName = 'SansSerif';
cb1.Label.FontSize = 11;

% b.2. nbar vs da & od 
nexttile
hm2 = imagesc(da_c, od_c, (Aggs_freq{2,2})');
xlim(del_da)
ylim(del_od)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('z_o_p_t (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb2 = colorbar;
cb2.Label.String = '$\bar{n}$';
cb2.Label.Interpreter = 'Latex';
cb2.Label.FontName = 'SansSerif';
cb2.Label.FontSize = 11;

% b.3. nbar vs ci & od
nexttile
hm3 = imagesc(ci_c, od_c, (Aggs_freq{2,3})');
xlim(del_ci)
ylim(del_od)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('z_o_p_t (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb3 = colorbar;
cb3.Label.String = '$\bar{n}$';
cb3.Label.Interpreter = 'Latex';
cb3.Label.FontName = 'SansSerif';
cb3.Label.FontSize = 11;

% c. the 3 parameter scatter plots

% c.1. od colormaps of ci vs da
nexttile
for k = 1 : n_bin(3)
    kk = Aggs_bin(:,3) == k;
    sp1 = scatter(da(kk), ci(kk), 10, cm3(k,:));
    hold on
end
axis equal
box on
xlim(del_da)
ylim(del_ci)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
legtxt1 = num2str(od_c, '%1.1f');
legend(legtxt1, 'Location', 'best', 'FontName', 'SansSerif',...
    'FontSize', 11)
sp1.Title.String = 'z_o_p_t (-)';
sp1.Title.FontName = 'SansSerif';
sp1.Title.FontSize = 11;

% c.2. ci colormaps of od vs da
nexttile
for j = 1 : n_bin(2)
    jj = Aggs_bin(:,2) == j;
    sp2 = scatter(da(jj), od(jj), 10, cm2(j,:));
    hold on
end
axis equal
box on
xlim(del_da)
ylim(del_od)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('z_o_p_t (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
legtxt2 = num2str(od_c, '%1.1f');
legend(legtxt2, 'Location', 'best', 'FontName', 'SansSerif',...
    'FontSize', 11)
sp2.Title.String = 'c (-)';
sp2.Title.FontName = 'SansSerif';
sp2.Title.FontSize = 11;

% c.3. da colormaps of od vs ci
nexttile
for i = 1 : n_bin(1)
    ii = Aggs_bin(:,1) == i;
    sp3 = scatter(ci(ii), od(kk), 10, cm1(i,:));
    hold on
end
axis equal
box on
xlim(del_ci)
ylim(del_od)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('d_a (m)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
legtxt3 = num2str(od_c, '%1.1f');
legend(legtxt3, 'Location', 'best', 'FontName', 'SansSerif',...
    'FontSize', 11)
sp3.Title.String = 'd_a (m)';
sp3.Title.FontName = 'SansSerif';
sp3.Title.FontSize = 11;

%%%


if nargout < 3
    clear h; % delete figure handle if not requested as an output
end

end

