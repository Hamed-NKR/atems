clc
clear
close all
warning('off')

%% initialize dpp vs. da figure %%

figure;
h = gcf;
h.Position = [100, 100, 500, 500];
set(h, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 2e1) ^ (1 / (1e4 - 1));
da0 = 2e1 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
dpp0 = 17.8 * (da0 / 100) .^ (0.35);
plt_0 = plot(da0, dpp0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

%% "low" agglomeration level data %%

% load files for the aggregate segmentation data
fname_agg_lal_1 = '19AUG-LAL-Start-Slider.mat';
fdir_agg_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ATEMS_Area';
id_agg_lal_1 = [5, 7, 10, 11, 27, 28, 30, 39, 42, 58, 60, 74];
fadd_agg_lal_1 = cell2mat(strcat(fdir_agg_lal_1, {'\'}, fname_agg_lal_1));
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
fname_pp_lal_1 = {'PFA_ET+NIT_HAL_20AUG24_R1082', 'PFA_ET+NIT_HAL_20AUG24_R1084',...
    'PFA_ET+NIT_HAL_20AUG24_R1087', 'PFA_ET+NIT_HAL_20AUG24_R1088',...
    'PFA_ET+NIT_HAL_20AUG24_R1107', 'PFA_ET+NIT_HAL_20AUG24_R1108',...
    'PFA_ET+NIT_HAL_20AUG24_R1110',...
    'PFA_ET+NIT_HAL_20AUG24_R1120', 'PFA_ET+NIT_HAL_20AUG24_R1123',...
    'PFA_ET+NIT_HAL_20AUG24_R1138', 'PFA_ET+NIT_HAL_20AUG24_R1140',...
    'PFA_ET+NIT_HAL_20AUG24_R1152'};
fdir_pp_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ImageJ_Primaries\19AUG24_LAL_R1\CSV';

n_agg_lal_1 = length(id_agg_lal_1);

dpp_manu_lal_1 = cell(n_agg_lal_1, 1);
dbarpp_manu_lal_1 = zeros(n_agg_lal_1, 1);
sigmapp_manu_lal_1 = zeros(n_agg_lal_1, 1);
npp_manu_lal_1 = zeros(n_agg_lal_1, 1);

fadd_pp_lal_1 = cell(n_agg_lal_1, 1);

for i = 1 : n_agg_lal_1

    fadd_pp_lal_1{i} = char(strcat(fdir_pp_lal_1, {'\'}, fname_pp_lal_1{i}, '.csv'));

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
        
    end

end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([20,2000])
ylim([10,60])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% "high" agglomeration level data %%

% load files for the aggregate segmentation data
fname_agg_hal_1 = '20AUG-HAL-Start-Slider.mat';
fdir_agg_hal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ATEMS_Area';
id_agg_hal_1 = [9, 11, 20, 48];
fadd_agg_hal_1 = cell2mat(strcat(fdir_agg_hal_1, {'\'}, fname_agg_hal_1));
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
fname_pp_hal_1 = {'20AUG24_HAL_9', '20AUG24_HAL_11', '20AUG24_HAL_20',...
    '20AUG24_HAL_48'};
fdir_pp_hal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ImageJ_Primaries\20AUG24_HAL_R1\CSV';

n_agg_hal_1 = length(id_agg_hal_1);

dpp_manu_hal_1 = cell(n_agg_hal_1, 1);
dbarpp_manu_hal_1 = zeros(n_agg_hal_1, 1);
sigmapp_manu_hal_1 = zeros(n_agg_hal_1, 1);
npp_manu_hal_1 = zeros(n_agg_hal_1, 1);

fadd_pp_hal_1 = cell(n_agg_hal_1, 1);

for i = 1 : n_agg_hal_1

    fadd_pp_hal_1{i} = char(strcat(fdir_pp_hal_1, {'\'}, fname_pp_hal_1{i}, '.csv'));

    if isfile(fadd_pp_hal_1{i})

        % opts_pp = detectImportOptions(fname_pp{i});
        % opts_pp = setvartype(opts_pp, 'char');
        pp_manu_hal_1 = readtable(fadd_pp_hal_1{i});
        
        dpp_manu_hal_1{i} = sqrt(4 * pp_manu_hal_1.Area(2:end) / pi);

        dbarpp_manu_hal_1(i) = geomean(dpp_manu_hal_1{i});
        sigmapp_manu_hal_1(i) = morph.geostd(dpp_manu_hal_1{i});
        npp_manu_hal_1(i) = length(dpp_manu_hal_1{i});
        
        plt_hal_1 = scatter(Aggs_hal_1(id_agg_hal_1(i)).da, dbarpp_manu_hal_1(i), 15,...
            hex2rgb('#C96868'), 's');
        
    end

end

legend(cat(2, plt_0, plt_lal_1, plt_hal_1), cat(2, {'Olfert and Rogak (2019)'},...
    {'Low agglomeration'}, {'High agglomeration'}),...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'northwest')

%% da comparison figure %%

figure;
h = gcf;
h.Position = [100, 100, 500, 500];
set(h, 'color', 'white');

condition = categorical([repmat({'Low agglomeration'}, size(Aggs_lal_1, 2), 1);...
    repmat({'High agglomeration'}, size(Aggs_hal_1, 2), 1)]);

boxplot([cat(1, Aggs_lal_1.da); cat(1, Aggs_hal_1.da)], condition)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
% ylim([10,60])
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

