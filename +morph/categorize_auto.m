function [fn, hf] = categorize_auto(Aggs, n_bin)

% CATEGORIZE_AUTO performs machine-based characterization of the...
%   ...particles based on their different morphological properties.
% ----------------------------------------------------------------------- %
% n_bin: numbers of post-processing bins (a 4 member vector for size,...
%   ...circularity, optical depth, and sharpness; also can be a single...
%   ...universal number for all)
% fn: A structure containing the morphological frequency data
% Aggs: aggregates' table of properties
% hf: A structure for figure handles of the frequency plots
% ----------------------------------------------------------------------- %

%% Initializations %%

n_agg_tot0 = length(Aggs); % total number of aggregates
fn0 = zeros(n_agg_tot0,5); % array to store bin labels
fn = struct(); % frequency storage structure
hf = struct(); % plot storage structure
    
% compile the aggregate (quantitative) properties
da = cat(1, Aggs.da); % area equivalent diameter
ca = cat(1, Aggs.ca); % area equivalent circularity
od = cat(1, Aggs.zbar_opt); % optical depth
% cat(1, Aggs.sbar_opt);
os = rand(n_agg_tot0,1); % optical sharpness

% reload the morphological types
for ind = 1 : n_agg_tot0
    if strcmp(Aggs(ind).Type, 'Fractal soot')
        fn0(ind,5) = 1;
    elseif strcmp(Aggs(ind).Type, 'Compact soot')
        fn0(ind,5) = 2;
    elseif strcmp(Aggs(ind).Type, 'Tarball')
        fn0(ind,5) = 3;
    elseif strcmp(Aggs(ind).Type, 'Softball')
        fn0(ind,5) = 4;
    elseif strcmp(Aggs(ind).Type, 'Hybrid')
        fn0(ind,5) = 5;
    elseif strcmp(Aggs(ind).Type, 'Miscellaneous')
        fn0(ind,5) = 6;
    else
        fn0(ind,5) = 0;
    end
end

% assign the number of bins if not given
if ~exist('n_bin', 'var') || isempty(n_bin)
    n_bin = [10, 10, 10, 10];
elseif length(n_bin) == 1
    n_bin = repmat(n_bin, 1, 4);
elseif length(n_bin(:)) ~= 4
    error(['Invalid bin nuumber set!', newline,...
        'Should be a 4 member vector'])
end

% range of properties
del_da = [min(da), max(da)];
del_ca = [min(ca), max(ca)];
del_od = [min(od), max(od)];
del_os = [min(os), max(os)];

% discretize the domains of post-processing variables (binning)
da_bin = del_da(1) + (del_da(2) - del_da(1)) .*...
    log(1 : (exp(1) - 1) / n_bin(1) : exp(1));
ca_bin = del_ca(1) + (del_ca(2) - del_ca(1)) .* (0 : 1 / n_bin(2) : 1);
od_bin = del_od(1) + (del_od(2) - del_od(1)) .* (0 : 1 / n_bin(3) : 1);
os_bin = del_os(1) + (del_os(2) - del_os(1)) .* (0 : 1 / n_bin(4) : 1);

% initialize the bin centers
da_c = zeros(n_bin(1),1);
ca_c = zeros(n_bin(2),1);
od_c = zeros(n_bin(3),1);
os_c = zeros(n_bin(4),1);

% Removing diverged/ data
ii0 = da > 1e4;
jj0 = (ca > 2) | (ca < -1);
kk0 = (od > 2) | (od < -1);
ll0 = (os > 2) | (os < -1);
mm0 = fn0(:,5) == 0;
rmv = ii0 | jj0 | kk0 | ll0 | mm0;
da(rmv) = [];
ca(rmv) = [];
od(rmv) = [];
os(rmv) = [];
fn0(rmv,:) = [];

n_agg_tot = length(da); % nubmer of aggregates after removing the...
    % ...diverged data

fn1 = {zeros(n_bin(1),1), zeros(n_bin(2),1), zeros(n_bin(3),1),...
    zeros(n_bin(4),1), zeros(6,1)}; % initialize the 1d number frequency...
        % ...array
fn2 = {zeros(n_bin(1), n_bin(2)), zeros(n_bin(1), n_bin(3)),...
    zeros(n_bin(1), n_bin(4)), zeros(n_bin(1), 6),...
    zeros(n_bin(2), n_bin(3)), zeros(n_bin(2), n_bin(4)),...
    zeros(n_bin(2), 6), zeros(n_bin(3), n_bin(4)),...
    zeros(n_bin(3), 6), zeros(n_bin(4), 6),}; % 2d frequencies



%% assign to 1d bins %%

for i = 1 : n_bin(1)
    if i == 1 % The first bin (only upper limit)
        ii = da < da_bin(i+1);
    elseif i < n_bin(1)
        ii = (da >= da_bin(i)) & (da < da_bin(i+1));
    else % the last bin (all the remaining)
        ii = fn0(:,1) == 0;
    end

    da_c(i) = sqrt(da_bin(i) * da_bin(i+1));
    fn0(ii,1) = i; % column 1 to be size bin labels
    fn1{1}(i) = nnz(ii) / n_agg_tot; % binned size freqs. (%)
end

% circularity binning
for j = 1 : n_bin(2)
    if j == 1
        jj = ca < ca_bin(j+1);
    elseif j < n_bin(2)
        jj = (ca >= ca_bin(j)) & (ca < ca_bin(j+1));
    else
        jj = fn0(:,2) == 0;
    end
    
    ca_c(j) = (ca_bin(j) + ca_bin(j+1)) / 2;
    fn0(jj,2) = j; % column 2 to be circularity bin labels    
    fn1{2}(j) = nnz(jj) / n_agg_tot; % circularity freqs. (%)
end

% optical depth binning
for k = 1 : n_bin(3)
    if k == 1
        kk = od < od_bin(k+1);
    elseif k < n_bin(3)
        kk = (od >= od_bin(k)) & (od < od_bin(k+1));
    else
        kk = fn0(:,3) == 0;
    end
    
    od_c(k) = (od_bin(k) + od_bin(k+1)) / 2;
    fn0(kk,3) = k; % column 3 to be optical depth bin labels
    fn1{3}(k) = nnz(kk) / n_agg_tot; % optical depth freqs. (%)
end

% optical sharpness binning
for l = 1 : n_bin(4)
    if l == 1
        ll = os < os_bin(l+1);
    elseif l < n_bin(4)
        ll = (os >= os_bin(l)) & (os < os_bin(l+1));
    else
        ll = fn0(:,4) == 0;
    end
    
    os_c(l) = (os_bin(l) + os_bin(l+1)) / 2;
    fn0(ll,4) = l; % column 4 to be optical sharpness bin labels
    fn1{4}(l) = nnz(ll) / n_agg_tot; % optical sharpness freqs. (%)
end

% morphological type binning
for m = 1 : 6
    fn1{5}(m) = nnz(fn0(:,5) == m) / n_agg_tot;
end

%% Plot univariate frequency barcharts %%

% initialize the figure
figure;
hf1 = gcf;
hf1.Position = [0, 0, 1500, 900]; % Position and size
set(hf1, 'color', 'white'); % Background color

% set the figure layout
tt1 = tiledlayout(2,3);
tt1.TileSpacing = 'loose';
tt1.Padding = 'loose';

% circularity subplot
nexttile(1);

bc1 = bar(ca_c, 100 * fn1{2}, 'FaceColor', 'flat');
cm11 = colormap(autumn);
jjj = round(1 + (length(cm11) - 1) .* (0.1 : 0.8 / (n_bin(2) - 1) : 0.9)');
cm11 = cm11(jjj,:); % get the descretized colormap
cm11 = flip(cm11,1);
bc1.CData = cm11;
hold on

err1 = sqrt((fn1{2} / n_agg_tot) .* (1 + fn1{2}) / n_agg_tot);
    % standard error propagation
eb1 = errorbar(ca_c, 100 * fn1{2}, 100 * err1, '.');
eb1.Color = [0.1 0.1 0.1];
eb1.CapSize = 5;
eb1.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.015 0.015], 'TickDir', 'out')

xticks(ca_c)
xtickformat('%.2f')
xtickangle(45)
ylim([0 10 * ceil(10 * max(fn1{2} + err1))])
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% optical depth subplot
nexttile(2);

bc2 = bar(od_c, 100 * fn1{3}, 'FaceColor', 'flat');
cm12 = colormap(gray);
kkk = round(1 + (length(cm12) - 1) .* (0.1 : 0.8 / (n_bin(3) - 1) : 0.9)');
cm12 = cm12(kkk,:);
cm12 = flip(cm12,1);
bc2.CData = cm12;
hold on

err2 = sqrt((fn1{3} / n_agg_tot) .* (1 + fn1{3}) / n_agg_tot);
eb2 = errorbar(od_c, 100 * fn1{3}, 100 * err2, '.');
eb2.Color = [0.6350 0.0780 0.1840];
eb2.CapSize = 5;
eb2.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.015 0.015], 'TickDir', 'out')
xticks(od_c)
xtickformat('%.2f')
xtickangle(45)
ylim([0 10 * ceil(10 * max(fn1{3} + err2))])
xlabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% sharpness subplot
nexttile(3);
bc3 = bar(os_c, 100 * fn1{4}, 'FaceColor', 'flat');
cm13 = colormap(winter);
lll = round(1 + (length(cm13) - 1) .* (0.1 : 0.8 / (n_bin(4) - 1) : 0.9)');
cm13 = cm13(lll,:);
bc3.CData = cm13;
hold on

err3 = sqrt((fn1{4} / n_agg_tot) .* (1 + fn1{4}) / n_agg_tot);
eb3 = errorbar(os_c, 100 * fn1{4}, 100 * err3, '.');
eb3.Color = [0.1 0.1 0.1];
eb3.CapSize = 5;
eb3.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.015 0.015], 'TickDir', 'out')
xticks(os_c)
xtickformat('%.2f')
xtickangle(45)
ylim([0 10 * ceil(10 * max(fn1{4} + err3))])
xlabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% morphological type subplot
nexttile(4);
x4 = categorical({'Fractal soot', 'Compact soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
x4 = reordercats(x4,{'Fractal soot', 'Compact soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
bc4 = bar(x4, 100 * fn1{5}, 'FaceColor', 'flat');
cm14 = colormap(summer);
mmm = round(1 + (length(cm14) - 1) .* (0.1 : 0.8 / 5 : 0.9)');
cm14 = cm14(mmm,:);
cm14 = flip(cm14,1);
bc4.CData = cm14;
hold on

err4 = sqrt((fn1{5} / n_agg_tot) .* (1 + fn1{5}) / n_agg_tot);
eb4 = errorbar(x4, 100 * fn1{5}, 100 * err4, '.');
eb4.Color = [0.1 0.1 0.1];
eb4.CapSize = 5;
eb4.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.015 0.015], 'TickDir', 'out')
xtickangle(45)
ylim([0 10 * ceil(10 * max(fn1{5} + err4))])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

title(tt1, 'Distributions of different morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
hold off

% size distribution subplot
nexttile(5, [1,2]);

fn11 = zeros(n_bin(1),1); % initialize dn/dlog(da)
err5 = zeros(n_bin(1),1); % initialize the errors
cm15 = colormap(spring);
iii = round(1 + (length(cm15) - 1) .* (0.1 : 0.8 / (n_bin(1) - 1) : 0.9)');
cm15 = cm15(iii,:);
for i = 1 : n_bin(1)
    fn11(i) = fn1{1}(i) / log(da_bin(i+1) / da_bin(i));
    err5(i) = sqrt((fn1{1}(i) / n_agg_tot) * (1 + fn1{1}(i))) /...
        (log(da_bin(i+1) / da_bin(i) * n_agg_tot));
    
    rectangle('Position',[da_bin(i), 0, (da_bin(i+1) - da_bin(i)),...
        fn11(i)], 'FaceColor', cm15(i,:));
    hold on
end
hold on

eb5 = errorbar(da_c, fn11, err5, '.');
eb5.Color = [0.1 0.1 0.1];
eb5.CapSize = 5;
eb5.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11, 'TickLength',...
    [0.0075 0.0075], 'TickDir', 'out')
set(gca, 'XScale', 'log')
set(gca, 'Layer', 'top')
xticks(da_bin)
xtickformat('%.0f')
xtickangle(45)
set(gca, 'XMinorTick', 'off')
xlim([del_da(1), del_da(2)])
ylim([0, ceil(10 * max(fn11 + err5)) / 10])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency size distribution (nm^-^1)',...
    'FontName', 'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
hold off

title(tt1, 'Univariate frequencies of morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

%% assign to 2d bins %%

% size distribution vs. circularity
for i = 1 : n_bin(1)
    for j = 1 : n_bin(2)
        fn2{1}(i,j) = nnz((fn0(:,1) == i) & (fn0(:,2) == j)) / n_agg_tot;
    end
end

% size distribution vs. optical depth
for i = 1 : n_bin(1)
    for k = 1 : n_bin(3)
        fn2{2}(i,k) = nnz((fn0(:,1) == i) & (fn0(:,3) == k)) / n_agg_tot;
    end
end

% size distribution vs. sharpness
for i = 1 : n_bin(1)
    for l = 1 : n_bin(4)
        fn2{3}(i,l) = nnz((fn0(:,1) == i) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% size distribution vs. morphological type
for i = 1 : n_bin(1)
    for m = 1 : 6
        fn2{4}(i,m) = nnz((fn0(:,1) == i) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% circularity vs. optical depth
for j = 1 : n_bin(2)
    for k = 1 : n_bin(3)
        fn2{5}(j,k) = nnz((fn0(:,2) == j) & (fn0(:,3) == k)) / n_agg_tot;
    end
end

% circularity vs. sharpness
for j = 1 : n_bin(2)
    for l = 1 : n_bin(4)
        fn2{6}(j,l) = nnz((fn0(:,2) == j) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% circularity vs. morphological type
for j = 1 : n_bin(2)
    for m = 1 : 6
        fn2{7}(j,m) = nnz((fn0(:,2) == j) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% optical depth vs. sharpness
for k = 1 : n_bin(3)
    for l = 1 : n_bin(4)
        fn2{8}(k,l) = nnz((fn0(:,3) == k) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% optical depth vs. morphological type
for k = 1 : n_bin(3)
    for m = 1 : 6
        fn2{9}(k,m) = nnz((fn0(:,3) == k) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% sharpness vs. morphological type
for l = 1 : n_bin(4)
    for m = 1 : 6
        fn2{10}(l,m) = nnz((fn0(:,4) == l) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

%% Plot bivariate frequency maps %%

% initialize the figure
figure;
hf2 = gcf;
hf2.Position = [0, 0, 1800, 900]; % Position and size
set(hf2, 'color', 'white'); % Background color

% set the figure layout
tt2 = tiledlayout(6,10);
tt2.TileSpacing = 'loose';
tt2.Padding = 'loose';

% size distribution vs circularity subplot
tt2_1 = nexttile(1, [2,2]);

pc1 = pcolor(da_bin, ca_bin, 100 * [(fn2{1})', zeros(n_bin(2), 1);...
    zeros(1, n_bin(1) + 1)]);
cm21 = colormap(tt2_1, hot);
cm21 = flip(cm21,1);
colormap(tt2_1, cm21)
% pc1.FaceColor = 'interp';
% pc1.EdgeColor = 'none';

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
yticks(ca_bin)
ytickformat('%.2f')
ylim(del_ca)
ylabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb1 = colorbar;
cb1.FontName = 'SansSerif';
cb1.FontSize = 10;
cb1.TickLength = 0.03;
cb1.Label.String = 'Frequency (%)';
cb1.Label.FontName = 'SansSerif';
cb1.Label.FontSize = 10;
cb1.Label.FontWeight = 'bold';
cb1.Label.Rotation = 270;
cb1_pos = get(cb1, 'Position');
cb1.Label.Position = [cb1_pos(1) + 3.5, cb1_pos(2) + 11]; % adjust...
    % ...colorbar label position

hold off

% size distribution vs optical depth subplot
tt2_2 = nexttile(21, [2,2]);

pc2 = pcolor(da_bin, od_bin, 100 * [(fn2{2})', zeros(n_bin(3), 1);...
    zeros(1, n_bin(1) + 1)]);
cm22 = colormap(tt2_2, gray);
cm22 = flip(cm22,1);
colormap(tt2_2, cm22)
% pc2.FaceColor = 'interp';
% pc2.EdgeColor = 'none';

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylim(del_od)
yticks(od_bin)
ytickformat('%.2f')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb2 = colorbar;
cb2.FontName = 'SansSerif';
cb2.FontSize = 10;
cb2.TickLength = 0.03;
cb2.Label.String = 'Frequency (%)';
cb2.Label.FontName = 'SansSerif';
cb2.Label.FontSize = 10;
cb2.Label.FontWeight = 'bold';
cb2.Label.Rotation = 270;
cb2_pos = get(cb2, 'Position');
cb2.Label.Position = [cb2_pos(1) + 3.5, cb2_pos(2) + 5.4];

hold off

% size distribution vs sharpness subplot
tt2_3 = nexttile(41, [2,2]);

pc3 = pcolor(da_bin, os_bin, 100 * [(fn2{3})', zeros(n_bin(4), 1);...
    zeros(1, n_bin(1) + 1)]);
colormap(tt2_3, parula)
% pc3.FaceColor = 'interp';
% pc3.EdgeColor = 'none';

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylim(del_os)
yticks(os_bin)
ytickformat('%.2f')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb3 = colorbar;
cb3.FontName = 'SansSerif';
cb3.FontSize = 10;
cb3.TickLength = 0.03;
cb3.Label.String = 'Frequency (%)';
cb3.Label.FontName = 'SansSerif';
cb3.Label.FontSize = 10;
cb3.Label.FontWeight = 'bold';
cb3.Label.Rotation = 270;
cb3_pos = get(cb3, 'Position');
cb3.Label.Position = [cb3_pos(1) + 3.5, cb3_pos(2) + 2.8];

hold off

% circularity vs optical depth subplot
tt2_4 = nexttile(3, [2,2]);

isc4 = imagesc(ca_c, od_c, 100 * (fn2{5})');
cm24 = colormap(tt2_4, copper);
cm24 = flip(cm24,1);
colormap(tt2_4, cm24)
% isc4.Interpolation = 'bilinear';

box on
axis padded
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1, 'YDir', 'normal')
xlim(del_ca)
xticks(ca_bin)
xtickformat('%.2f')
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylim(del_od)
yticks(od_bin)
ytickformat('%.2f')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb4 = colorbar;
cb4.FontName = 'SansSerif';
cb4.FontSize = 10;
cb4.TickLength = 0.03;
cb4.Label.String = 'Frequency (%)';
cb4.Label.FontName = 'SansSerif';
cb4.Label.FontSize = 10;
cb4.Label.FontWeight = 'bold';
cb4.Label.Rotation = 270;
cb4_pos = get(cb4, 'Position');
cb4.Label.Position = [cb4_pos(1) + 3.5, cb4_pos(2) + 5.2];

hold off

% circularity vs sharpness subplot
tt2_5 = nexttile(23, [2,2]);

isc5 = imagesc(ca_c, os_c, 100 * (fn2{6})');
cm25 = colormap(tt2_5, pink);
cm25 = flip(cm25,1);
colormap(tt2_5, cm25)
% isc5.Interpolation = 'bilinear';

box on
axis padded
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1, 'YDir', 'normal')
xticks(ca_bin)
xtickformat('%.2f')
xlim(del_ca)
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
yticks(os_bin)
ytickformat('%.2f')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb5 = colorbar;
cb5.FontName = 'SansSerif';
cb5.FontSize = 10;
cb5.TickLength = 0.03;
cb5.Label.String = 'Frequency (%)';
cb5.Label.FontName = 'SansSerif';
cb5.Label.FontSize = 10;
cb5.Label.FontWeight = 'bold';
cb5.Label.Rotation = 270;
cb5_pos = get(cb5, 'Position');
cb5.Label.Position = [cb5_pos(1) + 3.5, cb5_pos(2) + 3.3];

hold off

% optical depth vs sharpness subplot
tt2_6 = nexttile(43, [2,2]);

isc6 = imagesc(od_c, os_c, 100 * (fn2{8})');
cm26 = colormap(tt2_6, bone);
cm26 = flip(cm26,1);
colormap(tt2_6, cm26)
% isc6.Interpolation = 'bilinear';

box on
axis padded
grid on
xticks(od_bin)
xtickformat('%.2f')
xlim(del_od)
xlabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
yticks(os_bin)
ytickformat('%.2f')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1, 'YDir', 'normal')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb6 = colorbar;
cb6.FontName = 'SansSerif';
cb6.FontSize = 10;
cb6.TickLength = 0.03;
cb6.Label.String = 'Frequency (%)';
cb6.Label.FontName = 'SansSerif';
cb6.Label.FontSize = 10;
cb6.Label.FontWeight = 'bold';
cb6.Label.Rotation = 270;
cb6_pos = get(cb6, 'Position');
cb6.Label.Position = [cb6_pos(1) + 3.5, cb6_pos(2) + 2];

hold off

% size distribution vs morphlogical type
tt2_7 = nexttile(5, [3,3]);

x7 = categorical({'Fractal soot', 'Compact soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
x7 = reordercats(x7,{'Fractal soot', 'Compact soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
pc7 = pcolor((1 : 7), da_bin, 100 * [fn2{4}, zeros(n_bin(1), 1);...
    zeros(1, 7)]);
colormap(tt2_7, spring)
% pc7.FaceColor = 'interp';
% pc7.EdgeColor = 'none';

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1)
xticks(1.5 : 1 : 6.5)
xticklabels(x7)
xtickangle(45)
xlim([1,7])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
set(gca, 'YScale', 'log')
% yticks(da_bin)
% ytickformat('%.2f')
ylim(del_da)
ylabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
cb7 = colorbar;
cb7.FontName = 'SansSerif';
cb7.FontSize = 10;
cb7.TickLength = 0.015;
cb7.Label.String = 'Frequency (%)';
cb7.Label.FontName = 'SansSerif';
cb7.Label.FontSize = 10;
cb7.Label.FontWeight = 'bold';
cb7.Label.Rotation = 270;
cb7_pos = get(cb7, 'Position');
cb7.Label.Position = [cb7_pos(1) + 3.5, cb7_pos(2) + 3.6];

hold off

% circularity vs morphlogical type
tt2_8 = nexttile(35, [3,3]);

isc8 = pcolor((1 : 7), ca_bin, 100 * [fn2{7}, zeros(n_bin(2), 1);...
    zeros(1, 7)]);
cm28 = colormap(tt2_8, autumn);
cm28 = flip(cm28,1);
colormap(tt2_8, cm28)
% isc8.Interpolation = 'bilinear';

box on
axis padded
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1)
xticks(1.5 : 6.5)
xticklabels(x7)
xtickangle(45)
xlim([1, 7])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
yticks(ca_bin)
ytickformat('%.2f')
ylim(del_ca)
ylabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

cb8 = colorbar;
cb8.FontName = 'SansSerif';
cb8.FontSize = 10;
cb8.TickLength = 0.015;
cb8.Label.String = 'Frequency (%)';
cb8.Label.FontName = 'SansSerif';
cb8.Label.FontSize = 10;
cb8.Label.FontWeight = 'bold';
cb8.Label.Rotation = 270;
cb8_pos = get(cb8, 'Position');
cb8.Label.Position = [cb8_pos(1) + 3.5, cb8_pos(2) + 4.9];

hold off

% optical depth vs morphlogical type
tt2_9 = nexttile(8, [3,3]);

isc9 = pcolor((1 : 7), od_bin, 100 * [fn2{9}, zeros(n_bin(3), 1);...
    zeros(1, 7)]);
cm29 = colormap(tt2_9, gray);
cm29 = flip(cm29,1);
colormap(tt2_9, cm29)
% isc9.Interpolation = 'bilinear';

box on
axis padded
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1)
xticks(1.5 : 6.5)
xticklabels(x7)
xtickangle(45)
xlim([1, 7])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
yticks(od_bin)
ytickformat('%.2f')
ylim(del_od)
ylabel('Optical depth (-)', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
cb9 = colorbar;
cb9.FontName = 'SansSerif';
cb9.FontSize = 10;
cb9.TickLength = 0.015;
cb9.Label.String = 'Frequency (%)';
cb9.Label.FontName = 'SansSerif';
cb9.Label.FontSize = 10;
cb9.Label.FontWeight = 'bold';
cb9.Label.Rotation = 270;
cb9_pos = get(cb9, 'Position');
cb9.Label.Position = [cb9_pos(1) + 3.5, cb9_pos(2) + 3.6];

hold off

% optical sharpness vs morphlogical type
tt2_10 = nexttile(38, [3,3]);

isc10 = pcolor((1 : 7), os_bin, 100 * [fn2{10}, zeros(n_bin(4), 1);...
    zeros(1, 7)]);
colormap(tt2_10, winter)
% isc10.Interpolation = 'bilinear';

box on
axis padded
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out', 'GridColor', [0, 0, 0],...
    'GridAlpha', 1)
xticks(1.5 : 6.5)
xticklabels(x7)
xtickangle(45)
xlim([1, 7])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
yticks(os_bin)
ytickformat('%.2f')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')

cb10 = colorbar;
cb10.FontName = 'SansSerif';
cb10.FontSize = 10;
cb10.TickLength = 0.015;
cb10.Label.String = 'Frequency (%)';
cb10.Label.FontName = 'SansSerif';
cb10.Label.FontSize = 10;
cb10.Label.FontWeight = 'bold';
cb10.Label.Rotation = 270;
cb10_pos = get(cb10, 'Position');
cb10.Label.Position = [cb10_pos(1) + 3.5, cb10_pos(2) + 2.4];

hold off

title(tt2, 'Bivariate frequencies of morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

%% Draw trivariate colorcoded scatter plots %%

% initialize the figure
figure;
hf3 = gcf;
hf3.Position = [0, 0, 1500, 900]; % Position and size
set(hf3, 'color', 'white'); % Background color

% set the figure layout
tt3 = tiledlayout(2,3);
tt3.TileSpacing = 'loose';
tt3.Padding = 'loose';

% set the colormap
cm3 = colormap(turbo);
m4 = round(1 + (length(cm3) - 1) .* (0.1 : 0.8 / 5 : 0.9)');
cm3 = flip(cm3(m4,:), 1);
m5 = cell(6,1);

% size vs circularity
nexttile

for m = 1 : 6
    m5{m} = (fn0(:,5) == m);
    scatter(da(m5{m}), ca(m5{m}), 15, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylim(del_ca)
ylabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

% size vs optical depth
nexttile

for m = 1 : 6
    scatter(da(m5{m}), od(m5{m}), 15, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylim(del_od)
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

% size vs sharpness
nexttile

for m = 1 : 6
    scatter(da(m5{m}), os(m5{m}), 15, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
set(gca, 'XScale', 'log')
xlim(del_da)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

% circularity vs optical depth
nexttile

for m = 1 : 6
    scatter(ca(m5{m}), od(m5{m}), 15, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlim(del_ca)
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylim(del_od)
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

% circularity vs optical sharpness
nexttile

for m = 1 : 6
    scatter(ca(m5{m}), os(m5{m}), 10, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlim(del_ca)
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

% optical depth vs optical sharpness
nexttile

for m = 1 : 6
    scatter(od(m5{m}), os(m5{m}), 10, cm3(m,:))
    hold on
end

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlim(del_od)
xlabel('Optial depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylim(del_os)
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

hold off

legtxt3 = {'Fractal soot', 'Compact soot', 'Tarball', 'Softball',...
    'Hybrid', 'Miscellaneous'};
leg3 = legend(legtxt3, 'FontName', 'SansSerif', 'FontSize', 11);
leg3.Layout.Tile = 'east';
leg3.Title.String = 'Morphological type';
leg3.Title.FontName = 'SansSerif';
leg3.Title.FontSize = 12;
leg3.Title.FontWeight = 'bold';

title(tt3, 'Trivariate dispersions of morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

%% Store outputs %%

% update the frequency structure
fn.univar = fn1;
fn.bivar = fn2;
fn.labels = fn0;
fn.bins = {da_bin, ca_bin, od_bin, os_bin, x4};

% update the plot structure if figures are requested as output,...
    % ...otherwise delete
if nargout < 2
    clear hf1 hf2 hf3
else
    hf.univar = hf1;
    hf.bivar = hf2;
    hf.trivar = hf3;
end

end
