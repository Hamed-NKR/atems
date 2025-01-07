clc
clear
close all
warning('off')

%% initialize dpp vs. da figure %%

f1 = figure;
f1.Position = [50, 50, 600, 600];
set(f1, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 1e0) ^ (1 / (1e4 - 1));
da0 = 1e0 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
dpp0 = 17.8 * (da0 / 100) .^ (0.35);
plt_0 = plot(da0, dpp0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

%% "low" agglomeration level data - Repeat 1 %%

% load files for the aggregate segmentation data

fname_agg_lal_1 = '19AUG24_LAL_End_Slider';
fdir_agg_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ATEMS_Area';
fname_pp_lal_1 = 'PFA_ET+NIT_LAL_19AUG24_End';
fdir_pp_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ImageJ_Primaries\CSV';
id_agg_lal_1 = [1:27, 39, 42, 44, 55:64];

fadd_agg_lal_1 = cell2mat(strcat(fdir_agg_lal_1, {'\'}, fname_agg_lal_1, '.mat'));
load(fadd_agg_lal_1);

% Rename Aggs variable

vars = who; % call variable names

for i = 1 : length(vars)

    varname = vars{i};  % Get the name of the variable as a string

    if strcmp(varname, 'Aggs') % check if the variable contains Aggs
        newVarName = ['Aggs' '_lal_1']; % Modify the variable name by adding a suffix

        % Assign the value of the original variable to the new variable...
        %   ...(i.e. dynamically create a new variable)
        eval([newVarName ' = ' varname ';']);
        clear(varname); % Delete the old variable
    end

end

clear vars varname newVarName

% information on primary particle manual sizing data

n_agg_lal_1 = length(id_agg_lal_1);
dpp_manu_lal_1 = cell(n_agg_lal_1, 1);
dbarpp_manu_lal_1 = zeros(n_agg_lal_1, 1);
sigmapp_manu_lal_1 = zeros(n_agg_lal_1, 1);
npp_manu_lal_1 = zeros(n_agg_lal_1, 1);

fadd_pp_lal_1 = cell(n_agg_lal_1, 1);

figure(f1)
for i = 1 : n_agg_lal_1
    
    fadd_pp_lal_1{i} = char(strcat(fdir_pp_lal_1, '\', fname_pp_lal_1,...
        '_', num2str(int8(id_agg_lal_1(i))), '.csv'));

    if isfile(fadd_pp_lal_1{i})

        % opts_pp = detectImportOptions(fname_pp{i});
        % opts_pp = setvartype(opts_pp, 'char');
        pp_manu_lal_1 = readtable(fadd_pp_lal_1{i});
        
        dpp_manu_lal_1{i} = sqrt(4 * pp_manu_lal_1.Area(2:end) / pi);

        dbarpp_manu_lal_1(i) = geomean(dpp_manu_lal_1{i});
        sigmapp_manu_lal_1(i) = morph.geostd(dpp_manu_lal_1{i});
        npp_manu_lal_1(i) = length(dpp_manu_lal_1{i});
        
        plt_lal_1 = scatter(Aggs_lal_1(id_agg_lal_1(i)).da, dbarpp_manu_lal_1(i), 15,...
            hex2rgb('#006989'), '^');
        hold on
        
    end

end

%% "high" agglomeration level data - Repeat 1 %%

% load files for the aggregate segmentation data
fname_agg_hal_1 = '28AUG24_HAL_End_Slider';
fdir_agg_hal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\26SEP24_PFA_ET+NIT_HAL_28AUG24_End\ATEMS_Area';
fname_pp_hal_1 = 'PFA_ET+NIT_28AUG24_HAL_End';
fdir_pp_hal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\26SEP24_PFA_ET+NIT_HAL_28AUG24_End\ImageJ_Primaries\CSV';
id_agg_hal_1 = 1:40;

fadd_agg_hal_1 = cell2mat(strcat(fdir_agg_hal_1, {'\'}, fname_agg_hal_1, '.mat'));
load(fadd_agg_hal_1);

% Rename Aggs variable

vars = who; % call variable names

for i = 1 : length(vars)

    varname = vars{i};  % Get the name of the variable as a string

    if strcmp(varname, 'Aggs') % check if the variable contains Aggs
        newVarName = ['Aggs' '_hal_1']; % Modify the variable name by adding a suffix

        % Assign the value of the original variable to the new variable...
        %   ...(i.e. dynamically create a new variable)
        eval([newVarName ' = ' varname ';']);
        clear(varname); % Delete the old variable
    end

end

clear vars varname newVarName

% information on primary particle manual sizing data

n_agg_hal_1 = length(id_agg_hal_1);

dpp_manu_hal_1 = cell(n_agg_hal_1, 1);
dbarpp_manu_hal_1 = zeros(n_agg_hal_1, 1);
sigmapp_manu_hal_1 = zeros(n_agg_hal_1, 1);
npp_manu_hal_1 = zeros(n_agg_hal_1, 1);

fadd_pp_hal_1 = cell(n_agg_hal_1, 1);

for i = 1 : n_agg_hal_1

    fadd_pp_hal_1{i} = char(strcat(fdir_pp_hal_1, '\', fname_pp_hal_1,...
        '_', num2str(int8(id_agg_hal_1(i))), '.csv'));

    if isfile(fadd_pp_hal_1{i})

        % opts_pp = detectImportOptions(fname_pp{i});
        % opts_pp = setvartype(opts_pp, 'char');
        pp_manu_hal_1 = readtable(fadd_pp_hal_1{i});
        
        dpp_manu_hal_1{i} = sqrt(4 * pp_manu_hal_1.Area(2:end) / pi);

        dbarpp_manu_hal_1(i) = geomean(dpp_manu_hal_1{i});
        sigmapp_manu_hal_1(i) = morph.geostd(dpp_manu_hal_1{i});
        npp_manu_hal_1(i) = length(dpp_manu_hal_1{i});
        
        figure(f1);
        plt_hal_1 = scatter(Aggs_hal_1(id_agg_hal_1(i)).da, dbarpp_manu_hal_1(i), 15,...
            hex2rgb('#C96868'), 's');
        hold on
        
    end

end

% plot configs in dpp vs da figure
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([20,1000])
ylim([10,40])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

% number of aggregates to be manually sized
n_aggs_manu = [size(dpp_manu_lal_1, 1), size(dpp_manu_hal_1, 1)];

legend(cat(2, plt_lal_1, plt_hal_1, plt_0),...
    cat(2, strcat('Low agglom. (n =', {' '}, num2str(n_aggs_manu(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_aggs_manu(2)), ')'),...
    {'Olfert and Rogak (2019)'}), 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'northwest')


%% da comparison subplot %%

% initialize figure 2
f2 = figure;
f2.Position = [100, 100, 900, 900];
set(f2, 'color', 'white');

tt = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile

n_aggs_tot = [size(Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0), 2),...
    size(Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0), 2)];

xlbl21 = [strcat('Low agglom. (n =', {' '}, num2str(n_aggs_tot(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_aggs_tot(2)), ')')];

condition21 = [repmat(xlbl21(1), n_aggs_tot(1), 1);...
    repmat(xlbl21(2), n_aggs_tot(2), 1)];
condition21 = categorical(condition21, {xlbl21{1}, xlbl21{2}});

bp21 = boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)],...
    condition21, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.25);

boxes21 = findobj(bp21, 'Tag', 'Box');
patch(get(boxes21(1), 'XData'), get(boxes21(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes21(2), 'XData'), get(boxes21(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);

medians21 = findobj(bp21, 'Tag', 'Median');
set(medians21(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians21(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);

outliers21 = findobj(bp21, 'Tag', 'Outliers');
outliers21(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers21(1).MarkerSize = 3;
outliers21(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers21(2).MarkerSize = 3;

upwhisker21 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker21, 'linestyle', '-');
lowwhisker21= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker21, 'linestyle', '-');

hold on

% Compute kernel density estimate for each condition
[f_da_lal_1, xi_da_lal_1] =...
    ksdensity(log10(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da)));
[f_da_hal_1, xi_da_hal_1] =...
    ksdensity(log10(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)));

% to avoid issue with log scale in y axis
f_da_lal_1(xi_da_lal_1 <= 0) = [];
xi_da_lal_1(xi_da_lal_1 <= 0) = [];
f_da_hal_1(xi_da_hal_1 <= 0) = [];
xi_da_hal_1(xi_da_hal_1 <= 0) = [];

scale21 = -0.1;

plot(scale21 * f_da_lal_1 + 0.7, 10.^xi_da_lal_1, 'Color', hex2rgb('#8D493A'),...
    'LineWidth', 1.25)
fill([scale21 * f_da_lal_1 + 0.7, 0.7 * ones(size(f_da_lal_1))],...
     [10.^xi_da_lal_1, fliplr(10.^xi_da_lal_1)], hex2rgb('#DC8686'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

plot(scale21 * f_da_hal_1 + 1.7, 10.^xi_da_hal_1, 'Color', hex2rgb('#537188'),...
    'LineWidth', 1.25)
fill([scale21 * f_da_hal_1 + 1.7, 1.7 * ones(size(f_da_hal_1))],...
     [10.^xi_da_hal_1, fliplr(10.^xi_da_hal_1)], hex2rgb('#7EACB5'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
xlim([0.3, 2.7])
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
xlim([0.3, 2.7])
ylim([20, 1000])


%% ensemble dpp comparison subplot %%

nexttile

dpp_ens_lal_1 = cat(1, cell2mat(dpp_manu_lal_1));
dpp_ens_hal_1 = cat(1, cell2mat(dpp_manu_hal_1));

n_pps_manu = [size(dpp_ens_lal_1, 1), size(dpp_ens_hal_1, 1)];

xlbl22 = [strcat('Low agglom. (n =', {' '}, num2str(n_pps_manu(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_pps_manu(2)), ')')];

condition22 = [repmat(xlbl22(1), n_pps_manu(1), 1);...
    repmat(xlbl22(2), n_pps_manu(2), 1)];
condition22 = categorical(condition22, {xlbl22{1}, xlbl22{2}});

bp22 = boxplot([dpp_ens_lal_1; dpp_ens_hal_1], condition22, 'Notch','on',...
    'Symbol', 'o', 'Widths', 0.25);

boxes22 = findobj(bp22, 'Tag', 'Box');
patch(get(boxes22(1), 'XData'), get(boxes22(1), 'YData'),...
    hex2rgb('#DC8686'), 'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes22(2), 'XData'), get(boxes22(2), 'YData'),...
    hex2rgb('#7EACB5'), 'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);

medians22 = findobj(bp22, 'Tag', 'Median');
set(medians22(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians22(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);

outliers22 = findobj(bp22, 'Tag', 'Outliers');
outliers22(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers22(1).MarkerSize = 3;
outliers22(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers22(2).MarkerSize = 3;

upwhisker22 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker22, 'linestyle', '-');
lowwhisker22= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker22, 'linestyle', '-');

hold on

% Compute kernel density estimate for each condition
[f_dpp_lal_1, xi_dpp_lal_1] = ksdensity(log10(dpp_ens_lal_1));
[f_dpp_hal_1, xi_dpp_hal_1] = ksdensity(log10(dpp_ens_hal_1));

% to avoid issue with log scale in y axis
f_dpp_lal_1(xi_dpp_lal_1 <= 0) = [];
xi_dpp_lal_1(xi_dpp_lal_1 <= 0) = [];
f_dpp_hal_1(xi_dpp_hal_1 <= 0) = [];
xi_dpp_hal_1(xi_dpp_hal_1 <= 0) = [];

scale22 = -0.1;

plot(scale22 * f_dpp_lal_1 + 0.7, 10.^xi_dpp_lal_1, 'Color',...
    hex2rgb('#8D493A'), 'LineWidth', 1.25)
fill([scale22 * f_dpp_lal_1 + 0.7, 0.7 * ones(size(f_dpp_lal_1))],...
     [10.^xi_dpp_lal_1, fliplr(10.^xi_dpp_lal_1)], hex2rgb('#DC8686'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

plot(scale22 * f_dpp_hal_1 + 1.7, 10.^xi_dpp_hal_1, 'Color',...
    hex2rgb('#537188'), 'LineWidth', 1.25)
fill([scale22 * f_dpp_hal_1 + 1.7, 1.7 * ones(size(f_dpp_hal_1))],...
     [10.^xi_dpp_hal_1, fliplr(10.^xi_dpp_hal_1)], hex2rgb('#7EACB5'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
xlim([0.3, 2.7])
ylim([3, 80])

%% avereage dpp within aggregates comparison subplot %%

nexttile

xlbl23 = [strcat('Low agglom. (n =', {' '}, num2str(n_aggs_manu(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_aggs_manu(2)), ')')];

condition23_4 = categorical([repmat(xlbl23(1), n_aggs_manu(1), 1);...
    repmat(xlbl23(2), n_aggs_manu(2), 1)]);

bp23 = boxplot([dbarpp_manu_lal_1; dbarpp_manu_hal_1], condition23_4,...
    'Notch', 'on', 'Symbol', 'o', 'Widths', 0.3);

boxes23 = findobj(bp23, 'Tag', 'Box');
patch(get(boxes23(1), 'XData'), get(boxes23(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes23(2), 'XData'), get(boxes23(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);

medians23 = findobj(bp23, 'Tag', 'Median');
set(medians23(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians23(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);

outliers23 = findobj(bp23, 'Tag', 'Outliers');
outliers23(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers23(1).MarkerSize = 3;
outliers23(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers23(2).MarkerSize = 3;

upwhisker23 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker23, 'linestyle', '-');
lowwhisker23= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker23, 'linestyle', '-');

hold on

% ensemble geometric mean of primary particle size 
dbarpp_ens_lal_1 = geomean(dpp_ens_lal_1);
dbarpp_ens_hal_1 = geomean(dpp_ens_hal_1);
plot([1.6, 2.4], [dbarpp_ens_lal_1, dbarpp_ens_lal_1], 'Color',...
    [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':')
plt23_ens = plot([0.6, 1.4], [dbarpp_ens_hal_1, dbarpp_ens_hal_1], 'Color',...
    [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');

legend(plt23_ens, '$\langle{d_\mathrm{pp}}\rangle$ (Ensemble GM)',...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southeast')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'XDir', 'reverse')
ylabel('$\overline{d}_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylim([11.5, 29.5])

%% GSD of pp within aggregates comparison subplot %%

nexttile

bp24 = boxplot([sigmapp_manu_lal_1; sigmapp_manu_hal_1], condition23_4,...
    'Notch', 'on', 'Symbol', 'o', 'Widths', 0.3);
hold on

boxes24 = findobj(bp24, 'Tag', 'Box');
patch(get(boxes24(1), 'XData'), get(boxes24(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes24(2), 'XData'), get(boxes24(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);

medians24 = findobj(bp24, 'Tag', 'Median');
set(medians24(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians24(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);

outliers24 = findobj(bp24, 'Tag', 'Outliers');
outliers24(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers24(1).MarkerSize = 3;
outliers24(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers24(2).MarkerSize = 3;

upwhisker24 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker24, 'linestyle', '-');
lowwhisker24= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker24, 'linestyle', '-');

hold on

% ensemble geometric standard deviation of primary particle size 
sigma_ens_lal_1 = morph.geostd(dpp_ens_lal_1);
sigma_ens_hal_1 = morph.geostd(dpp_ens_hal_1);

plot([1.6, 2.4], [sigma_ens_lal_1, sigma_ens_lal_1], 'Color',...
    [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':')
plt24_ens = plot([0.6, 1.4], [sigma_ens_hal_1, sigma_ens_hal_1], 'Color',...
    [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');

legend(plt24_ens, '$\gamma_\mathrm{pp}$ (Ensemble GSD)', 'interpreter',...
    'latex', 'FontSize', 11, 'location', 'southeast')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log', 'XDir', 'reverse')
yticks([1.2 1.3 1.4 1.5 1.6])
ylabel('$\sigma_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylim([1.15, 1.65])

%% n_hyb freqeuncy comparison subplot %%

% initialize figure 3
f3 = figure;
f3.Position = [150, 150, 500, 400];
set(f3, 'color', 'white');

n_subagg_0 = {cat(1,Aggs_lal_1.n_subagg), cat(1,Aggs_hal_1.n_subagg)};

n_subagg = n_subagg_0;

n_subagg{1}(n_subagg_0{1} < 1) = [];
n_subagg{2}(n_subagg_0{2} < 1) = [];

n_hyb = 100 * [nnz(n_subagg{1} == 1),...
    nnz(n_subagg{1} == 2),...
    nnz((n_subagg{1} >= 3) & (n_subagg{1} <= 5)),...
    nnz(n_subagg{1} > 5);...
    nnz(n_subagg{2} == 1),...
    nnz(n_subagg{2} == 2),...
    nnz((n_subagg{2} >= 3) & (n_subagg{2} <= 5)),...
    nnz(n_subagg{2} > 5)];
n_hyb(1,:) = n_hyb(1,:) / length(n_subagg{1});
n_hyb(2,:) = n_hyb(2,:) / length(n_subagg{2});
    
condition3 = categorical(xlbl21);

b3 = bar(condition3, n_hyb,'stacked');
clr3 = hex2rgb({'#4793AF', '#FFC470', '#DD5746', '#8B322C'});

for i = 1 : size(n_hyb,2)
    b3(i).BarWidth = 0.2;
    b3(i).FaceColor = clr3(i,:);
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02])
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('Frequency [$\%$]', 'interpreter', 'latex', 'FontSize', 14)

lgd3 = legend({'$n_\mathrm{hyb} = 1$', '$n_\mathrm{hyb} = 2$',...
    '$3 \le n_\mathrm{hyb} \le 5$', '$n_\mathrm{hyb} > 5$'},...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'northoutside',...
    'orientation', 'horizontal', 'NumColumns', 4);
lgd3.ItemTokenSize = [10, 10];

%% dpp vs. da as a function of n_hyb, not experimental case %%

% initialize figure 4
f4 = figure;
f4.Position = [200, 200, 600, 600];
set(f4, 'color', 'white')

plt4 = cell(5,1); % initialize dpp vs da plots per n_hyb

 % compile number of subaggregates
n_subagg_tot = cat(1, n_subagg_0{1}(id_agg_lal_1), n_subagg_0{2}(id_agg_hal_1));

% compile mean primary particle sizes
dbarpp_tot = [dbarpp_manu_lal_1; dbarpp_manu_hal_1];

% compile projected-area sizes
da_tot = [cat(1, Aggs_lal_1(id_agg_lal_1).da);...
    cat(1, Aggs_hal_1(id_agg_hal_1).da)];

plt4{1} = plot(da0, dpp0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

plt4{2} = scatter(da_tot(n_subagg_tot == 1), dbarpp_tot(n_subagg_tot == 1),...
    15, clr3(1,:), '^');

plt4{3} = scatter(da_tot(n_subagg_tot == 2), dbarpp_tot(n_subagg_tot == 2),...
    15, clr3(2,:), 's');

plt4{4} = scatter(da_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)),...
    dbarpp_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)), 25,...
    clr3(3,:), 'h');

plt4{5} = scatter(da_tot(n_subagg_tot > 5), dbarpp_tot(n_subagg_tot > 5),...
    15, clr3(4,:), 'o');

% plot configs
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([20,1000])
ylim([10,40])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

% number of aggregates to be manually sized
n_aggs_manu = [size(dpp_manu_lal_1, 1), size(dpp_manu_hal_1, 1)];

legend(cat(2, plt4{:}),...
    {'Olfert and Rogak (2019)', '$n_\mathrm{hyb} = 1$',...
    '$n_\mathrm{hyb} = 2$', '$3 \le n_\mathrm{hyb} \le 5$',...
    '$n_\mathrm{hyb} > 5$'}, 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'northwest')

%% outputs for Langevin dynamics modeling

% geometric mean and standard deviation of aggregate projected area size...
    % ...for uniform-looking low-agglomeration aggregates
da_uni = cat(1, Aggs_lal_1.da);
da_uni(cat(1,Aggs_lal_1.n_subagg) ~= 1) = [];
GM_da_uni = geomean(da_uni);
GSD_da_uni = morph.geostd(da_uni);

% find indices of aggregate with uniform structure
ii0_agg_uni = find(cat(1,Aggs_lal_1.n_subagg) == 1);
[ii_agg_uni, ii_pp_uni] = intersect(id_agg_lal_1, ii0_agg_uni);

% geometric mean and standard deviation of geometric mean of primary...
    % ...particle diamter within individual aggregates (GM and GSD of mean)
dpp_uni = dbarpp_manu_lal_1(ii_pp_uni);
GM_dpp_uni = geomean(dpp_uni);
GSD_dpp_uni = morph.geostd(dpp_uni);

