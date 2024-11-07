clc
clear
close all
warning('off')

%% initialize dpp vs. da figure %%

f1 = figure;
f1.Position = [50, 50, 500, 600];
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
id_agg_hal_1 = 1:30;

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
    'location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2)


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

condition21 = categorical([repmat(xlbl21(1), n_aggs_tot(1), 1);...
    repmat(xlbl21(2), n_aggs_tot(2), 1)]);

boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)],...
    condition21, 'Notch', 'on')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% ensemble dpp comparison subplot %%

nexttile

dpp_ens_lal_1 = cat(1, cell2mat(dpp_manu_lal_1));
dpp_ens_hal_1 = cat(1, cell2mat(dpp_manu_hal_1));

n_pps_manu = [size(dpp_ens_lal_1, 1), size(dpp_ens_hal_1, 1)];

xlbl22 = [strcat('Low agglom. (n =', {' '}, num2str(n_pps_manu(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_pps_manu(2)), ')')];

condition22 = categorical([repmat(xlbl22(1), n_pps_manu(1), 1);...
    repmat(xlbl22(2), n_pps_manu(2), 1)]);

boxplot([dpp_ens_lal_1; dpp_ens_hal_1], condition22, 'Notch','on')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% avereage dpp within aggregates comparison subplot %%

nexttile

xlbl23 = [strcat('Low agglom. (n =', {' '}, num2str(n_aggs_manu(1)), ')'),...
    strcat('High agglom. (n =', {' '}, num2str(n_aggs_manu(2)), ')')];

condition23_4 = categorical([repmat(xlbl23(1), n_aggs_manu(1), 1);...
    repmat(xlbl23(2), n_aggs_manu(2), 1)]);

boxplot([dbarpp_manu_lal_1; dbarpp_manu_hal_1], condition23_4, 'Notch', 'on')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\overline{d}_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% GSD of pp within aggregates comparison subplot %%

nexttile

boxplot([sigmapp_manu_lal_1; sigmapp_manu_hal_1], condition23_4, 'Notch', 'on')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\sigma_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)

%% n_hyb freqeuncy comparison subplot %%

% initialize figure 3
f3 = figure;
f3.Position = [150, 150, 500, 400];
set(f3, 'color', 'white');

n_subagg = {cat(1,Aggs_lal_1.n_subagg), cat(1,Aggs_hal_1.n_subagg)};
n_subagg{1}(n_subagg{1} < 1) = [];
n_subagg{2}(n_subagg{2} < 1) = [];

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
f4.Position = [200, 200, 1200, 900];
set(f4, 'color', 'white')

plt4 = cell(5,1); % initialize dpp vs da plots per n_hyb

n_subagg_tot = cat(1, n_subagg{:}); % compile number of subaggregates
dbarpp_tot = [dbarpp_manu_lal_1; dbarpp_manu_hal_1]; % compile mean primary particle sizes
dbarpp_tot(n_subagg{1} < 1) = [];
da_tot = [cat(1, Aggs_lal_1.da); cat(1, Aggs_hal_1.da)]; % compile projected-area sizes

plt4{1} = copyobj(plt_0, f4.CurrentAxes);

plt4{2} = scatter(da_tot(n_subagg_tot == 1), dbarpp_tot(n_subagg_tot == 1),...
    15, hex2rgb(clr3(1,:)), '^');
hold on

plt4{3} = scatter(da_tot(n_subagg_tot == 2), dbarpp_tot(n_subagg_tot == 2),...
    15, hex2rgb(clr3(2,:)), 's');

plt4{4} = scatter(da_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)),...
    dbarpp_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)), 15,...
    hex2rgb(clr3(3,:)), 'h');

plt4{5} = scatter(da_tot(n_subagg_tot > 5), dbarpp_tot(n_subagg_tot > 5),...
    15, hex2rgb(clr3(4,:)), 'o');

% plot configs
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
    'location', 'eastoutside')
