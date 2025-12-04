clc
clear
close all
warning('off')

%% universal parameters %%


n_rndsmp = 5000; % number of random points for resampling needed to...
    % ...generate boxplots of ensemble primary particle distributions

% Logistic switch parameter
T = 1.1;   % threshold coverage
k = 0.05;   % bandwidth of transition

%% initialize dpp vs. da figure %%

f1 = figure;
f1.Position = [50, 50, 650, 800];
set(f1, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 1e0) ^ (1 / (1e4 - 1));
da0 = 1e0 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
dpp0 = 17.8 * (da0 / 100) .^ (0.35);
plt_0 = plot(da0, dpp0, 'Color', hex2rgb('#DEAA79'),... % [0.4940 0.1840 0.5560]
    'LineStyle', '-.', 'LineWidth', 3);
hold on

%% "low" agglomeration level data - Repeat 1 %%

% load files for the aggregate segmentation data

fname_agg_lal_1 = '19AUG24_LAL_End_Slider';
fdir_agg_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ATEMS_Area';
fname_pp_lal_1 = 'PFA_ET+NIT_LAL_19AUG24_End';
fdir_pp_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ImageJ_Primaries\CSV';
id_agg_lal_1 = 1:64;

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

        % coverage ratio = (sum PP areas) / (agg area) ; constants cancel on d^2
        da_i_lal = Aggs_lal_1(id_agg_lal_1(i)).da;
        cov_i_lal = sum(dpp_manu_lal_1{i}.^2) / (da_i_lal^2);
        
        % if cov_i_lal < T
        %     w_i_lal = 1 ./ (dpp_manu_lal_1{i}.^2); % apply correction
        % else
        %     w_i_lal = ones(size(dpp_manu_lal_1{i})); % no weighting
        % end
        alpha_i = 1 - 1 ./ (1 + exp(-(cov_i_lal - T)/k));
        w_i_lal = (1 ./ (dpp_manu_lal_1{i}.^2)) .^ alpha_i;

        [dbarpp_manu_lal_1(i), sigmapp_manu_lal_1(i)] = ...
            morph.geomstats(dpp_manu_lal_1{i}, w_i_lal);
        
        npp_manu_lal_1(i) = length(dpp_manu_lal_1{i});
        
        plt_lal_1 = scatter(Aggs_lal_1(id_agg_lal_1(i)).da,...
            dbarpp_manu_lal_1(i), 25, hex2rgb('#C96868'), '^',...
            'LineWidth', 1.5);
        hold on
        
    end

end

%% Load data of "Moderate" undiluted agglomeration condition - Aggregates have moderate hybridity and high collapse %%

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

        da_i_hal  = Aggs_hal_1(id_agg_hal_1(i)).da;
        cov_i_hal = sum(dpp_manu_hal_1{i}.^2) / (da_i_hal^2);
        % w_i_hal = (cov_i_hal < T) .* (1./(dpp_manu_hal_1{i}.^2)) +...
        %     (cov_i_hal >= T) .* 1;
        alpha_i = 1 - 1 ./ (1 + exp(-(cov_i_hal - T)/k));
        w_i_hal = (1 ./ (dpp_manu_hal_1{i}.^2)) .^ alpha_i;        
        [dbarpp_manu_hal_1(i), sigmapp_manu_hal_1(i)] =...
            morph.geomstats(dpp_manu_hal_1{i}, w_i_hal);
        
        npp_manu_hal_1(i) = length(dpp_manu_hal_1{i});
        
        figure(f1);
        plt_hal_1 = scatter(Aggs_hal_1(id_agg_hal_1(i)).da,...
            dbarpp_manu_hal_1(i), 25, hex2rgb('#006989'), 'o',...
            'LineWidth', 1.5);
        hold on
        
    end

end


%% Load data of high agglomeration diluted condition - Aggregates have high hybridity and low collapse %%

% load files for the aggregate segmentation data

fname_agg_exdil = '06FEB25_ExAglom_Slider';
fdir_agg_exdil = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\Extreme\ExAglom\Aggs\Data';
fname_pp_exdil = 'PFA_ET+NIT_06FEB25_EXAGLOM';
fdir_pp_exdil= 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\Extreme\ExAglom\PP\CSV';
id_agg_exdil = [1,5,9,10,12:18,21,22,25,29];

fadd_agg_exdil = cell2mat(strcat(fdir_agg_exdil, {'\'}, fname_agg_exdil, '.mat'));
load(fadd_agg_exdil);

% Rename Aggs variable

vars = who; % call variable names

for i = 1 : length(vars)

    varname = vars{i};  % Get the name of the variable as a string

    if strcmp(varname, 'Aggs') % check if the variable contains Aggs
        newVarName = ['Aggs' '_exdil']; % Modify the variable name by adding a suffix

        % Assign the value of the original variable to the new variable...
        %   ...(i.e. dynamically create a new variable)
        eval([newVarName ' = ' varname ';']);
        clear(varname); % Delete the old variable
    end

end

clear vars varname newVarName

% load primary particle size data
n_agg_exdil = length(id_agg_exdil);
dpp_manu_exdil = cell(n_agg_exdil, 1);
dbarpp_manu_exdil = zeros(n_agg_exdil, 1);
sigmapp_manu_exdil = zeros(n_agg_exdil, 1);
npp_manu_exdil = zeros(n_agg_exdil, 1);

fadd_pp_exdil = cell(n_agg_exdil, 1);

for i = 1 : n_agg_exdil

    fadd_pp_exdil{i} = char(strcat(fdir_pp_exdil, '\', fname_pp_exdil,...
        '_', num2str(int8(id_agg_exdil(i))), '.csv'));

    if isfile(fadd_pp_exdil{i})

        % opts_pp = detectImportOptions(fname_pp{i});
        % opts_pp = setvartype(opts_pp, 'char');
        pp_manu_exdil = readtable(fadd_pp_exdil{i});
        
        dpp_manu_exdil{i} = sqrt(4 * pp_manu_exdil.Area(2:end) / pi);

        da_i_exdil  = Aggs_exdil(id_agg_exdil(i)).da;
        cov_i_exdil = sum(dpp_manu_exdil{i}.^2) / (da_i_exdil^2);
        % w_i_exdil = (cov_i_exdil < T) .* (1./(dpp_manu_exdil{i}.^2)) +...
        %     (cov_i_exdil >= T) .* 1;
        alpha_i = 1 - 1 ./ (1 + exp(-(cov_i_exdil - T)/k));
        w_i_exdil = (1 ./ (dpp_manu_exdil{i}.^2)) .^ alpha_i;
        [dbarpp_manu_exdil(i), sigmapp_manu_exdil(i)] =...
            morph.geomstats(dpp_manu_exdil{i}, w_i_exdil);

        npp_manu_exdil(i) = length(dpp_manu_exdil{i});
        
        figure(f1);
        plt_exdil = scatter(Aggs_exdil(id_agg_exdil(i)).da,...
            dbarpp_manu_exdil(i), 35, hex2rgb('#8174A0'), 's',...
            'LineWidth', 1.5);
        hold on
        
    end

end

% number of aggregates to be manually sized
n_aggs_manu = [size(dpp_manu_lal_1, 1), size(dpp_manu_hal_1, 1),...
    size(dpp_manu_exdil, 1)];

% plot configs in dpp vs da figure
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
[xmin1, xmax1] = bounds(cat(1, cat(1, Aggs_lal_1.da), cat(1, Aggs_hal_1.da)));
[ymin1, ymax1] = bounds(cat(1, dbarpp_manu_lal_1, dbarpp_manu_hal_1));
xlim([0.85 * xmin1, 1.2 * xmax1])
ylim([0.95 * ymin1, 1.05 * ymax1])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
xticks(cat(2,linspace(30,100,8),linspace(200,700,6)))
xtickangle(90)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
legend(cat(2, plt_lal_1, plt_hal_1, plt_exdil, plt_0),...
    cat(2, strcat('Lo-Aglom (n =', {' '}, num2str(n_aggs_manu(1)), ')'),...
    strcat('Mod-Colaps (n =', {' '}, num2str(n_aggs_manu(2)), ')'),...
    strcat('Ex-Aglom (n =', {' '}, num2str(n_aggs_manu(3)), ')'),...
    {'Olfert and Rogak (2019)'}), 'interpreter', 'latex', 'FontSize', 16,...
    'location', 'southoutside', 'Orientation', 'horizontal',...
    'NumColumns', 2)

%% Load data of extreme undiluted agglomeration condition - Aggregates have high hybridity and high collapse %%

% load files for the aggregate segmentation data

fname_agg_excol = '07FEB25_ExColaps_Slider';
fdir_agg_excol = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\Extreme\ExColaps\Aggs\Data';

fadd_agg_excol = cell2mat(strcat(fdir_agg_excol, {'\'}, fname_agg_excol, '.mat'));
load(fadd_agg_excol);

% Rename Aggs variable

vars = who; % call variable names

for i = 1 : length(vars)

    varname = vars{i};  % Get the name of the variable as a string

    if strcmp(varname, 'Aggs') % check if the variable contains Aggs
        newVarName = ['Aggs' '_excol']; % Modify the variable name by adding a suffix

        % Assign the value of the original variable to the new variable...
        %   ...(i.e. dynamically create a new variable)
        eval([newVarName ' = ' varname ';']);
        clear(varname); % Delete the old variable
    end

end

clear vars varname newVarName


%% ensemble dpp comparison subplot %%

% initialize figure 2
f2 = figure;
f2.Position = [100, 0, 850, 1300];
set(f2, 'color', 'white');

tt2 = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile

% strcat('Low agglom.', string(newline), '(n =', {' '}, num2str(n_aggs_tot(1)), ')'),...
%     strcat('High agglom.', string(newline), '(n =', {' '}, num2str(n_aggs_tot(2)), ')')];
% 

% compile primary particle diamaters across aggregates
dpp_ens_lal_1 = cat(1, cell2mat(dpp_manu_lal_1));
dpp_ens_hal_1 = cat(1, cell2mat(dpp_manu_hal_1));
dpp_ens_exdil = cat(1, cell2mat(dpp_manu_exdil));

%%% apply weights to resolve bias in sizing of larger primary particles

% coverage per aggregate (PP area / agg area)
cov_ens_lal = nan(n_agg_lal_1,1);
for ii = 1 : n_agg_lal_1
    if ~isempty(dpp_manu_lal_1{ii})
        da_temp_ens = Aggs_lal_1(id_agg_lal_1(ii)).da;
        cov_ens_lal(ii) = sum(dpp_manu_lal_1{ii}.^2) / (da_temp_ens^2);
    end
end

cov_ens_hal = nan(n_agg_hal_1,1);
for ii = 1 : n_agg_hal_1
    if ~isempty(dpp_manu_hal_1{ii})
        da_temp_ens = Aggs_hal_1(id_agg_hal_1(ii)).da;
        cov_ens_hal(ii) = sum(dpp_manu_hal_1{ii}.^2) / (da_temp_ens^2);
    end
end

cov_ens_exdil = nan(n_agg_exdil,1);
for ii = 1 : n_agg_exdil
    if ~isempty(dpp_manu_exdil{ii})
        da_temp_ens = Aggs_exdil(id_agg_exdil(ii)).da;
        cov_ens_exdil(ii) = sum(dpp_manu_exdil{ii}.^2) / (da_temp_ens^2);
    end
end

% expand parent agg diameters and coverage to per-primary particle vectors
da_parent_lal = repelem(cat(1, Aggs_lal_1(id_agg_lal_1).da), npp_manu_lal_1, 1);
covvec_lal = repelem(cov_ens_lal, npp_manu_lal_1, 1);
da_parent_hal = repelem(cat(1, Aggs_hal_1(id_agg_hal_1).da), npp_manu_hal_1, 1);
covvec_hal = repelem(cov_ens_hal, npp_manu_hal_1, 1);
da_parent_exdil = repelem(cat(1, Aggs_exdil(id_agg_exdil).da), npp_manu_exdil, 1);
covvec_exdil = repelem(cov_ens_exdil, npp_manu_exdil, 1);

% ensemble weights: (A_agg / A_pp) = (da_parent.^2 ./ dpp.^2)

% w_lal = ones(size(dpp_ens_lal_1));
% mask = covvec_lal < T;
% w_lal(mask) = (da_parent_lal(mask).^2) ./ (dpp_ens_lal_1(mask).^2);
% 
% w_hal = ones(size(dpp_ens_hal_1));
% mask = covvec_hal < T;
% w_hal(mask) = (da_parent_hal(mask).^2) ./ (dpp_ens_hal_1(mask).^2);
% 
% w_exdil = ones(size(dpp_ens_exdil));
% mask = covvec_exdil < T;
% w_exdil(mask) = (da_parent_exdil(mask).^2) ./ (dpp_ens_exdil(mask).^2);

alpha_lal = 1 - 1 ./ (1 + exp(-(covvec_lal - T)/k));
w_lal = ((da_parent_lal.^2) ./ (dpp_ens_lal_1.^2)) .^ alpha_lal;

alpha_hal = 1 - 1 ./ (1 + exp(-(covvec_hal - T)/k));
w_hal = ((da_parent_hal.^2) ./ (dpp_ens_hal_1.^2)) .^ alpha_hal;

alpha_exdil = 1 - 1 ./ (1 + exp(-(covvec_exdil - T)/k));
w_exdil = ((da_parent_exdil.^2) ./ (dpp_ens_exdil.^2)) .^ alpha_exdil;

% weighted ensemble GM & GSD
[gm_lal, gsd_lal] = morph.geomstats(dpp_ens_lal_1, w_lal);
[gm_hal, gsd_hal] = morph.geomstats(dpp_ens_hal_1, w_hal);
[gm_exdil, gsd_exdil] = morph.geomstats(dpp_ens_exdil, w_exdil);

%%%

% generate weights to correct ensemble particle sizing bias

n_pps_manu = [size(dpp_ens_lal_1, 1), size(dpp_ens_hal_1, 1),...
    size(dpp_ens_exdil, 1)];

xlbl21 = [strcat('Lo-Aglom (n =', {' '}, num2str(n_pps_manu(1)), ')'),...
    strcat('Mod-Colaps (n =', {' '}, num2str(n_pps_manu(2)), ')'),...
    strcat('Hi-Aglom (n =', {' '}, num2str(n_pps_manu(3)), ')')];

condition21 = [repmat(xlbl21(1), n_pps_manu(1), 1);...
    repmat(xlbl21(2), n_pps_manu(2), 1);...
    repmat(xlbl21(3), n_pps_manu(3), 1)];
condition21 = categorical(condition21, {xlbl21{1}, xlbl21{2}, xlbl21{3}});

% resample data with weights and pass to boxplot
dpp_all = [dpp_ens_lal_1; dpp_ens_hal_1; dpp_ens_exdil];
w_all = [w_lal(:); w_hal(:); w_exdil(:)];
idx = randsample(numel(dpp_all), n_rndsmp, true, w_all/sum(w_all));
bp21 = boxplot(dpp_all(idx), condition21(idx), 'Notch', 'on',...
     'Symbol', 'o', 'Widths', 0.25);

boxes21 = findobj(bp21, 'Tag', 'Box');
patch(get(boxes21(1), 'XData'), get(boxes21(1), 'YData'),...
    hex2rgb('#DC8686'), 'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes21(2), 'XData'), get(boxes21(2), 'YData'),...
    hex2rgb('#7EACB5'), 'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes21(3), 'XData'), get(boxes21(3), 'YData'),...
    hex2rgb('#A888B5'), 'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);

medians21 = findobj(bp21, 'Tag', 'Median');
set(medians21(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians21(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians21(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);

outliers21 = findobj(bp21, 'Tag', 'Outliers');
outliers21(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers21(1).MarkerSize = 3;
outliers21(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers21(2).MarkerSize = 3;
outliers21(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers21(3).MarkerSize = 3;

upwhisker21 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker21, 'linestyle', '-');
lowwhisker21= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker21, 'linestyle', '-');

hold on

% Compute kernel density estimate for each condition
[f_dpp_lal_1, xi_dpp_lal_1] = ksdensity(log10(dpp_ens_lal_1),...
    'Weights', w_lal/sum(w_lal));
[f_dpp_hal_1, xi_dpp_hal_1] = ksdensity(log10(dpp_ens_hal_1),...
    'Weights', w_hal/sum(w_hal));
[f_dpp_exdil, xi_dpp_exdil] = ksdensity(log10(dpp_ens_exdil),...
    'Weights', w_exdil/sum(w_exdil));

% to avoid issue with log scale in y axis
f_dpp_lal_1(xi_dpp_lal_1 <= 0) = [];
xi_dpp_lal_1(xi_dpp_lal_1 <= 0) = [];
f_dpp_hal_1(xi_dpp_hal_1 <= 0) = [];
xi_dpp_hal_1(xi_dpp_hal_1 <= 0) = [];
f_dpp_exdil(xi_dpp_exdil <= 0) = [];
xi_dpp_exdil(xi_dpp_exdil <= 0) = [];

scale21 = -0.1;

plot(scale21 * f_dpp_lal_1 + 0.7, 10.^xi_dpp_lal_1, 'Color',...
    hex2rgb('#8D493A'), 'LineWidth', 1.25)
fill([scale21 * f_dpp_lal_1 + 0.7, 0.7 * ones(size(f_dpp_lal_1))],...
     [10.^xi_dpp_lal_1, fliplr(10.^xi_dpp_lal_1)], hex2rgb('#DC8686'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

plot(scale21 * f_dpp_hal_1 + 1.7, 10.^xi_dpp_hal_1, 'Color',...
    hex2rgb('#537188'), 'LineWidth', 1.25)
fill([scale21 * f_dpp_hal_1 + 1.7, 1.7 * ones(size(f_dpp_hal_1))],...
     [10.^xi_dpp_hal_1, fliplr(10.^xi_dpp_hal_1)], hex2rgb('#7EACB5'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

plot(scale21 * f_dpp_exdil + 2.7, 10.^xi_dpp_exdil, 'Color',...
    hex2rgb('#8174A0'), 'LineWidth', 1.25)
fill([scale21 * f_dpp_exdil + 2.7, 2.7 * ones(size(f_dpp_exdil))],...
     [10.^xi_dpp_exdil, fliplr(10.^xi_dpp_exdil)], hex2rgb('#A888B5'),...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.015 0.015], 'YScale', 'log')
ylabel('$d_\mathrm{pp}^\mathrm{(i)}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 24)
xlim([0.3, 3.3])
ylim([5, 55])
yticks([5 10 20 40 80])

%% avereage dpp within aggregates comparison subplot %%

nexttile

xlbl22 = [strcat('Lo-Aglom (n =', {' '}, num2str(n_aggs_manu(1)), ')'),...
    strcat('Mod-Colaps (n =', {' '}, num2str(n_aggs_manu(2)), ')'),...
    strcat('Hi-Aglom (n =', {' '}, num2str(n_aggs_manu(3)), ')')];

condition22_23 = categorical([repmat(xlbl22(1), n_aggs_manu(1), 1);...
    repmat(xlbl22(2), n_aggs_manu(2), 1);...
    repmat(xlbl22(3), n_aggs_manu(3), 1)]);
condition22_23 = categorical(condition22_23, {xlbl22{1}, xlbl22{2},...
    xlbl22{3}});

bp22 = boxplot([dbarpp_manu_lal_1; dbarpp_manu_hal_1; dbarpp_manu_exdil],...
    condition22_23, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.3);

boxes22 = findobj(bp22, 'Tag', 'Box');
patch(get(boxes22(1), 'XData'), get(boxes22(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes22(2), 'XData'), get(boxes22(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes22(3), 'XData'), get(boxes22(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);

medians22 = findobj(bp22, 'Tag', 'Median');
set(medians22(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians22(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians22(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);

outliers22 = findobj(bp22, 'Tag', 'Outliers');
outliers22(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers22(1).MarkerSize = 3;
outliers22(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers22(2).MarkerSize = 3;
outliers22(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers22(3).MarkerSize = 3;

upwhisker22 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker22, 'linestyle', '-');
lowwhisker22= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker22, 'linestyle', '-');

hold on

% ensemble geometric mean of primary particle size 
dbarpp_ens_lal_1 = geomean(dpp_ens_lal_1);
dbarpp_ens_hal_1 = geomean(dpp_ens_hal_1);
dbarpp_ens_exdil = geomean(dpp_ens_exdil);
% plot([0.6, 1.4], [dbarpp_ens_lal_1, dbarpp_ens_lal_1], 'Color',...
%     [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% plot([1.6, 2.4], [dbarpp_ens_hal_1, dbarpp_ens_hal_1], 'Color',...
%     [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% plt22_ens = plot([2.6, 3.4], [dbarpp_ens_exdil, dbarpp_ens_exdil],...
%     'Color', [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% 
% legend(plt22_ens, '$\langle{d_\mathrm{pp}}\rangle$ (Ensemble GM)',...
%     'interpreter', 'latex', 'FontSize', 10, 'location', 'northeast')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.015 0.015])
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
ylim([8, 24])

%% GSD of pp within aggregates comparison subplot %%

nexttile

bp23 = boxplot([sigmapp_manu_lal_1; sigmapp_manu_hal_1; sigmapp_manu_exdil],...
    condition22_23, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.3);
hold on

boxes23 = findobj(bp23, 'Tag', 'Box');
patch(get(boxes23(1), 'XData'), get(boxes23(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes23(2), 'XData'), get(boxes23(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes23(3), 'XData'), get(boxes23(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);

medians23 = findobj(bp23, 'Tag', 'Median');
set(medians23(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians23(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians23(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);

outliers23 = findobj(bp23, 'Tag', 'Outliers');
outliers23(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers23(1).MarkerSize = 3;
outliers23(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers23(2).MarkerSize = 3;
outliers23(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers23(3).MarkerSize = 3;

upwhisker23 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker23, 'linestyle', '-');
lowwhisker23= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker23, 'linestyle', '-');

hold on

% ensemble geometric standard deviation of primary particle size 
sigma_ens_lal_1 = morph.geostd(dpp_ens_lal_1);
sigma_ens_hal_1 = morph.geostd(dpp_ens_hal_1);
sigma_ens_exdil = morph.geostd(dpp_ens_exdil);

% plot([0.6, 1.4], [sigma_ens_lal_1, sigma_ens_lal_1], 'Color',...
%     [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% plot([1.6, 2.4], [sigma_ens_hal_1, sigma_ens_hal_1], 'Color',...
%     [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% plt23_ens = plot([2.6, 3.4], [sigma_ens_exdil, sigma_ens_exdil], 'Color',...
%     [0, 0, 0], 'LineWidth', 1.25, 'LineStyle', ':');
% 
% legend(plt23_ens, '$\gamma_\mathrm{pp}$ (Ensemble GSD)', 'interpreter',...
%     'latex', 'FontSize', 10, 'location', 'northeast')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.015 0.015])
yticks([1.2 1.3 1.4 1.5 1.6])
ylabel('$\sigma_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 24)
ylim([1.11, 1.62])

%% Summary table for Figure 2 (Geometric mean, GSD & 95% CI) %%
% all CIs are log-space normal approximations.

%% Summary table for Figure 2 (Ensemble variants by rows + aggregate summaries)
% all CIs are log-space normal approximations.

% ---------- 1) Build the three weighting modes for ENSEMBLE d_pp^(i) ----------
% Unweighted
w_lal_unw   = ones(size(dpp_ens_lal_1));
w_hal_unw   = ones(size(dpp_ens_hal_1));
w_exdil_unw = ones(size(dpp_ens_exdil));

% Hard threshold (A_agg/A_pp if coverage < T, else 1)
w_lal_hard   = ones(size(dpp_ens_lal_1));
mask         = covvec_lal < T;
w_lal_hard(mask) = (da_parent_lal(mask).^2) ./ (dpp_ens_lal_1(mask).^2);

w_hal_hard   = ones(size(dpp_ens_hal_1));
mask         = covvec_hal < T;
w_hal_hard(mask) = (da_parent_hal(mask).^2) ./ (dpp_ens_hal_1(mask).^2);

w_exdil_hard = ones(size(dpp_ens_exdil));
mask         = covvec_exdil < T;
w_exdil_hard(mask) = (da_parent_exdil(mask).^2) ./ (dpp_ens_exdil(mask).^2);

% Soft (logistic) are already computed earlier as: w_lal, w_hal, w_exdil
% If not present, uncomment these lines:
% alpha_lal   = 1 - 1 ./ (1 + exp(-(covvec_lal   - T)/k));
% alpha_hal   = 1 - 1 ./ (1 + exp(-(covvec_hal   - T)/k));
% alpha_exdil = 1 - 1 ./ (1 + exp(-(covvec_exdil - T)/k));
% w_lal   = ((da_parent_lal  .^2) ./ (dpp_ens_lal_1.^2)) .^ alpha_lal;
% w_hal   = ((da_parent_hal  .^2) ./ (dpp_ens_hal_1.^2)) .^ alpha_hal;
% w_exdil = ((da_parent_exdil.^2) ./ (dpp_ens_exdil.^2)) .^ alpha_exdil;

% ---------- 2) Compute ensemble stats for each method & condition ----------
% Unweighted
[gm_ens_lal_unw,  gsd_ens_lal_unw,  ci_ens_lal_unw]  = morph.geomstats(dpp_ens_lal_1(:),  w_lal_unw(:));
[gm_ens_hal_unw,  gsd_ens_hal_unw,  ci_ens_hal_unw]  = morph.geomstats(dpp_ens_hal_1(:),  w_hal_unw(:));
[gm_ens_exdil_unw,gsd_ens_exdil_unw,ci_ens_exdil_unw]= morph.geomstats(dpp_ens_exdil(:),  w_exdil_unw(:));

% Hard
[gm_ens_lal_hard,  gsd_ens_lal_hard,  ci_ens_lal_hard]  = morph.geomstats(dpp_ens_lal_1(:),  w_lal_hard(:));
[gm_ens_hal_hard,  gsd_ens_hal_hard,  ci_ens_hal_hard]  = morph.geomstats(dpp_ens_hal_1(:),  w_hal_hard(:));
[gm_ens_exdil_hard,gsd_ens_exdil_hard,ci_ens_exdil_hard]= morph.geomstats(dpp_ens_exdil(:),  w_exdil_hard(:));

% Soft (your current logistic results already computed as gm_* / gsd_* / ci_*)
[gm_ens_lal_soft,  gsd_ens_lal_soft,  ci_ens_lal_soft]  = morph.geomstats(dpp_ens_lal_1(:),  w_lal(:));
[gm_ens_hal_soft,  gsd_ens_hal_soft,  ci_ens_hal_soft]  = morph.geomstats(dpp_ens_hal_1(:),  w_hal(:));
[gm_ens_exdil_soft,gsd_ens_exdil_soft,ci_ens_exdil_soft]= morph.geomstats(dpp_ens_exdil(:),  w_exdil(:));

% ---------- 3) Assemble long-format ENSEMBLE table (rows = condition × method) ----------
Cond = categorical(repmat({'Lo-Aglom','Mod-Colaps','Hi-Aglom'}, 1, 3))';
Meth = categorical([repmat({'Unweighted'},3,1); repmat({'Hard'},3,1); repmat({'Soft'},3,1)]);

GM_ens  = [gm_ens_lal_unw;    gm_ens_hal_unw;    gm_ens_exdil_unw; ...
           gm_ens_lal_hard;   gm_ens_hal_hard;   gm_ens_exdil_hard; ...
           gm_ens_lal_soft;   gm_ens_hal_soft;   gm_ens_exdil_soft];

GSD_ens = [gsd_ens_lal_unw;   gsd_ens_hal_unw;   gsd_ens_exdil_unw; ...
           gsd_ens_lal_hard;  gsd_ens_hal_hard;  gsd_ens_exdil_hard; ...
           gsd_ens_lal_soft;  gsd_ens_hal_soft;  gsd_ens_exdil_soft];

CI_lo_ens = [ci_ens_lal_unw(1);  ci_ens_hal_unw(1);  ci_ens_exdil_unw(1); ...
             ci_ens_lal_hard(1); ci_ens_hal_hard(1); ci_ens_exdil_hard(1); ...
             ci_ens_lal_soft(1); ci_ens_hal_soft(1); ci_ens_exdil_soft(1)];

CI_hi_ens = [ci_ens_lal_unw(2);  ci_ens_hal_unw(2);  ci_ens_exdil_unw(2); ...
             ci_ens_lal_hard(2); ci_ens_hal_hard(2); ci_ens_exdil_hard(2); ...
             ci_ens_lal_soft(2); ci_ens_hal_soft(2); ci_ens_exdil_soft(2)];

T_fig2_ens_long = table(Cond, Meth, GM_ens, GSD_ens, CI_lo_ens, CI_hi_ens, ...
    'VariableNames', {'Condition','Method','GM_dpp_Ens','GSD_dpp_Ens','CI95_lo_dpp_Ens','CI95_hi_dpp_Ens'});

% ---------- 4) Keep your per-aggregate blocks as-is (unweighted across aggs) ----------
% subplot 2: per-aggregate mean d_pp
[gm_gmagg_lal, gsd_gmagg_lal, ci_gmagg_lal] = morph.geomstats(dbarpp_manu_lal_1);
[gm_gmagg_hal, gsd_gmagg_hal, ci_gmagg_hal] = morph.geomstats(dbarpp_manu_hal_1);
[gm_gmagg_exdil, gsd_gmagg_exdil, ci_gmagg_exdil] = morph.geomstats(dbarpp_manu_exdil);

T_fig2_gmagg = table( ...
    categorical({'Lo-Aglom';'Mod-Colaps';'Hi-Aglom'}), ...
    [gm_gmagg_lal;  gm_gmagg_hal;  gm_gmagg_exdil], ...
    [gsd_gmagg_lal; gsd_gmagg_hal; gsd_gmagg_exdil], ...
    [ci_gmagg_lal(1); ci_gmagg_hal(1); ci_gmagg_exdil(1)], ...
    [ci_gmagg_lal(2); ci_gmagg_hal(2); ci_gmagg_exdil(2)], ...
    'VariableNames', {'Condition','GM_dpp_GM','GSD_dpp_GM','CI95_lo_dpp_GM','CI95_hi_dpp_GM'});

% subplot 3: per-aggregate sigma_pp
[gm_gsdagg_lal, gsd_gsdagg_lal, ci_gsdagg_lal] = morph.geomstats(sigmapp_manu_lal_1);
[gm_gsdagg_hal, gsd_gsdagg_hal, ci_gsdagg_hal] = morph.geomstats(sigmapp_manu_hal_1);
[gm_gsdagg_exdil, gsd_gsdagg_exdil, ci_gsdagg_exdil] = morph.geomstats(sigmapp_manu_exdil);

T_fig2_gsdagg = table( ...
    categorical({'Lo-Aglom';'Mod-Colaps';'Hi-Aglom'}), ...
    [gm_gsdagg_lal;  gm_gsdagg_hal;  gm_gsdagg_exdil], ...
    [gsd_gsdagg_lal; gsd_gsdagg_hal; gsd_gsdagg_exdil], ...
    [ci_gsdagg_lal(1); ci_gsdagg_hal(1); ci_gsdagg_exdil(1)], ...
    [ci_gsdagg_lal(2); ci_gsdagg_hal(2); ci_gsdagg_exdil(2)], ...
    'VariableNames', {'Condition','GM_dpp_GSD','GSD_dpp_GSD','CI95_lo_dpp_GSD','CI95_hi_dpp_GSD'});

% ---------- 5) (Optional) show all three tables nicely ----------
disp(' ');
disp('=== Figure 2 — Ensemble d_pp^(i) by Method (rows) ===');
disp(T_fig2_ens_long);

disp(' ');
disp('=== Figure 2 — Per-aggregate mean d_pp (unweighted across aggregates) ===');
disp(T_fig2_gmagg);

disp(' ');
disp('=== Figure 2 — Per-aggregate sigma_pp (unweighted across aggregates) ===');
disp(T_fig2_gsdagg);

%% da comparison subplot %%

% initialize figure 3
f3 = figure;
f3.Position = [150, 50, 900, 900];
set(f3, 'color', 'white');

tt3 = tiledlayout(2, 2, 'Padding', 'loose', 'TileSpacing', 'compact');
nexttile

n_aggs_tot = [size(Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0), 2),...
    size(Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0), 2),...
    size(Aggs_exdil, 2), size(Aggs_excol, 2)];

xlbl31 = [strcat('Lo-Aglom', {' '}, '(n =', {' '}, num2str(n_aggs_tot(1)), ')'),...
    strcat('Mod-Colaps', {' '},  '(n =', {' '}, num2str(n_aggs_tot(2)), ')'),...
    strcat('Hi-Aglom', {' '},  '(n =', {' '}, num2str(n_aggs_tot(3)), ')'),...
    strcat('Ex-Colaps', {' '},  '(n =', {' '}, num2str(n_aggs_tot(4)), ')')];

% strcat('Low agglom.', string(newline), '(n =', {' '}, num2str(n_aggs_tot(1)), ')'),...
%     strcat('High agglom.', string(newline), '(n =', {' '}, num2str(n_aggs_tot(2)), ')')];
% 
condition31 = [repmat(xlbl31(1), n_aggs_tot(1), 1);...
    repmat(xlbl31(2), n_aggs_tot(2), 1);...
    repmat(xlbl31(3), n_aggs_tot(3), 1);...
    repmat(xlbl31(4), n_aggs_tot(4), 1)];
condition31 = categorical(condition31, {xlbl31{1}, xlbl31{2}, xlbl31{3},...
    xlbl31{4}});

bp31 = boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da);...
    cat(1, Aggs_exdil.da); cat(1, Aggs_excol.da)],...
    condition31, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.25);

GMs_da = [geomean(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da)),...
    geomean(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)),...
    geomean(cat(1, Aggs_exdil.da)), geomean(cat(1, Aggs_excol.da))];

GSDs_da = [morph.geostd(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da)),...
    morph.geostd(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)),...
    morph.geostd(cat(1, Aggs_exdil.da)), morph.geostd(cat(1, Aggs_excol.da))];

boxes31 = findobj(bp31, 'Tag', 'Box');
patch(get(boxes31(1), 'XData'), get(boxes31(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes31(2), 'XData'), get(boxes31(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes31(3), 'XData'), get(boxes31(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);
patch(get(boxes31(4), 'XData'), get(boxes31(4), 'YData'), hex2rgb('#B1C29E'),...
    'EdgeColor', hex2rgb('#659287'), 'FaceAlpha', 0.5);

medians31 = findobj(bp31, 'Tag', 'Median');
set(medians31(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians31(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians31(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);
set(medians31(4), 'Color', hex2rgb('#5F6F65'), 'LineWidth', 2);

outliers31 = findobj(bp31, 'Tag', 'Outliers');
outliers31(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers31(1).MarkerSize = 3;
outliers31(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers31(2).MarkerSize = 3;
outliers31(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers31(3).MarkerSize = 3;
outliers31(4).MarkerEdgeColor = hex2rgb('#B1C29E');
outliers31(4).MarkerSize = 3;

upwhisker31 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker31, 'linestyle', '-');
lowwhisker31= findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker31, 'linestyle', '-');

% hold on
% 
% % Compute kernel density estimate for each condition
% [f_da_lal_1, xi_da_lal_1] =...
%     ksdensity(log10(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).da)));
% [f_da_hal_1, xi_da_hal_1] =...
%     ksdensity(log10(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).da)));
% [f_da_exdil, xi_da_exdil] =...
%     ksdensity(log10(cat(1, Aggs_exdil.da)));
% [f_da_excol, xi_da_excol] =...
%     ksdensity(log10(cat(1, Aggs_excol.da)));
% 
% % to avoid issue with log scale in y axis
% f_da_lal_1(xi_da_lal_1 <= 0) = [];
% xi_da_lal_1(xi_da_lal_1 <= 0) = [];
% f_da_hal_1(xi_da_hal_1 <= 0) = [];
% xi_da_hal_1(xi_da_hal_1 <= 0) = [];
% f_da_exdil(xi_da_exdil <= 0) = [];
% xi_da_exdil(xi_da_exdil <= 0) = [];
% f_da_excol(xi_da_excol <= 0) = [];
% xi_da_excol(xi_da_excol <= 0) = [];
% 
% scale31 = -0.15;
% 
% plot(scale31 * f_da_lal_1 + 0.7, 10.^xi_da_lal_1, 'Color', hex2rgb('#8D493A'),...
%     'LineWidth', 1.25)
% fill([scale31 * f_da_lal_1 + 0.7, 0.7 * ones(size(f_da_lal_1))],...
%      [10.^xi_da_lal_1, fliplr(10.^xi_da_lal_1)], hex2rgb('#DC8686'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale31 * f_da_hal_1 + 1.7, 10.^xi_da_hal_1, 'Color', hex2rgb('#537188'),...
%     'LineWidth', 1.25)
% fill([scale31 * f_da_hal_1 + 1.7, 1.7 * ones(size(f_da_hal_1))],...
%      [10.^xi_da_hal_1, fliplr(10.^xi_da_hal_1)], hex2rgb('#7EACB5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale31 * f_da_exdil + 2.7, 10.^xi_da_exdil, 'Color', hex2rgb('#8174A0'),...
%     'LineWidth', 1.25)
% fill([scale31 * f_da_exdil + 2.7, 2.7 * ones(size(f_da_exdil))],...
%      [10.^xi_da_exdil, fliplr(10.^xi_da_exdil)], hex2rgb('#A888B5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale31 * f_da_excol + 3.7, 10.^xi_da_excol, 'Color', hex2rgb('#659287'),...
%     'LineWidth', 1.25)
% fill([scale31 * f_da_excol + 3.7, 3.7 * ones(size(f_da_excol))],...
%      [10.^xi_da_excol, fliplr(10.^xi_da_excol)], hex2rgb('#B1C29E'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
% xlim([0.3, 4.3])
ylim([20, 3000])

%% circularity distribution

nexttile

bp32 = boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).ca);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).ca);...
    cat(1, Aggs_exdil.ca); cat(1, Aggs_excol.ca)],...
    condition31, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.25);

boxes32 = findobj(bp32, 'Tag', 'Box');
patch(get(boxes32(1), 'XData'), get(boxes32(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes32(2), 'XData'), get(boxes32(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes32(3), 'XData'), get(boxes32(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);
patch(get(boxes32(4), 'XData'), get(boxes32(4), 'YData'), hex2rgb('#B1C29E'),...
    'EdgeColor', hex2rgb('#659287'), 'FaceAlpha', 0.5);

medians32 = findobj(bp32, 'Tag', 'Median');
set(medians32(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians32(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians32(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);
set(medians32(4), 'Color', hex2rgb('#5F6F65'), 'LineWidth', 2);

outliers32 = findobj(bp32, 'Tag', 'Outliers');
outliers32(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers32(1).MarkerSize = 3;
outliers32(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers32(2).MarkerSize = 3;
outliers32(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers32(3).MarkerSize = 3;
outliers32(4).MarkerEdgeColor = hex2rgb('#B1C29E');
outliers32(4).MarkerSize = 3;

upwhisker32 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker32, 'linestyle', '-');
lowwhisker32 = findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker32, 'linestyle', '-');

% hold on
% 
% % Compute kernel density estimate for each condition
% [f_ca_lal_1, xi_ca_lal_1] =...
%     ksdensity(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).ca));
% [f_ca_hal_1, xi_ca_hal_1] =...
%     ksdensity(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).ca));
% [f_ca_exdil, xi_ca_exdil] =...
%     ksdensity(cat(1, Aggs_exdil.ca));
% [f_ca_excol, xi_ca_excol] =...
%     ksdensity(cat(1, Aggs_excol.ca));
% 
% % to avoid issue with log scale in y axis
% f_ca_lal_1(xi_ca_lal_1 <= 0) = [];
% xi_ca_lal_1(xi_ca_lal_1 <= 0) = [];
% f_ca_hal_1(xi_ca_hal_1 <= 0) = [];
% xi_ca_hal_1(xi_ca_hal_1 <= 0) = [];
% f_ca_exdil(xi_ca_exdil <= 0) = [];
% xi_ca_exdil(xi_ca_exdil <= 0) = [];
% f_ca_excol(xi_ca_excol <= 0) = [];
% xi_ca_excol(xi_ca_excol <= 0) = [];
% 
% scale32 = -0.075;
% 
% plot(scale32 * f_ca_lal_1 + 0.7, xi_ca_lal_1, 'Color', hex2rgb('#8D493A'),...
%     'LineWidth', 1.25)
% fill([scale32 * f_ca_lal_1 + 0.7, 0.7 * ones(size(f_ca_lal_1))],...
%      [xi_ca_lal_1, fliplr(xi_ca_lal_1)], hex2rgb('#DC8686'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale32 * f_ca_hal_1 + 1.7, xi_ca_hal_1, 'Color', hex2rgb('#537188'),...
%     'LineWidth', 1.25)
% fill([scale32 * f_ca_hal_1 + 1.7, 1.7 * ones(size(f_ca_hal_1))],...
%      [xi_ca_hal_1, fliplr(xi_ca_hal_1)], hex2rgb('#7EACB5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale32 * f_ca_exdil + 2.7, xi_ca_exdil, 'Color', hex2rgb('#8174A0'),...
%     'LineWidth', 1.25)
% fill([scale32 * f_ca_exdil + 2.7, 2.7 * ones(size(f_ca_exdil))],...
%      [xi_ca_exdil, fliplr(xi_ca_exdil)], hex2rgb('#A888B5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale32 * f_ca_excol + 3.7, xi_ca_excol, 'Color', hex2rgb('#659287'),...
%     'LineWidth', 1.25)
% fill([scale32 * f_ca_excol + 3.7, 3.7 * ones(size(f_ca_excol))],...
%      [xi_ca_excol, fliplr(xi_ca_excol)], hex2rgb('#B1C29E'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.02 0.02])
ylabel('$c_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
% xlim([0.3, 4.3])
ylim([0, 1])

%% optical depth distribution

nexttile

bp33 = boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).zbar_opt);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).zbar_opt);...
    cat(1, Aggs_exdil.zbar_opt); cat(1, Aggs_excol.zbar_opt)],...
    condition31, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.25);

boxes33 = findobj(bp33, 'Tag', 'Box');
patch(get(boxes33(1), 'XData'), get(boxes33(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes33(2), 'XData'), get(boxes33(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes33(3), 'XData'), get(boxes33(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);
patch(get(boxes33(4), 'XData'), get(boxes33(4), 'YData'), hex2rgb('#B1C29E'),...
    'EdgeColor', hex2rgb('#659287'), 'FaceAlpha', 0.5);

medians33 = findobj(bp33, 'Tag', 'Median');
set(medians33(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians33(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians33(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);
set(medians33(4), 'Color', hex2rgb('#5F6F65'), 'LineWidth', 2);

outliers33 = findobj(bp33, 'Tag', 'Outliers');
outliers33(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers33(1).MarkerSize = 3;
outliers33(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers33(2).MarkerSize = 3;
outliers33(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers33(3).MarkerSize = 3;
outliers33(4).MarkerEdgeColor = hex2rgb('#B1C29E');
outliers33(4).MarkerSize = 3;

upwhisker33 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker33, 'linestyle', '-');
lowwhisker33 = findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker33, 'linestyle', '-');

% hold on
% 
% % Compute kernel density estimate for each condition
% [f_za_lal_1, xi_za_lal_1] =...
%     ksdensity(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).zbar_opt));
% [f_za_hal_1, xi_za_hal_1] =...
%     ksdensity(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).zbar_opt));
% [f_za_exdil, xi_za_exdil] =...
%     ksdensity(cat(1, Aggs_exdil.zbar_opt));
% [f_za_excol, xi_za_excol] =...
%     ksdensity(cat(1, Aggs_excol.zbar_opt));
% 
% % to avoid issue with log scale in y axis
% f_za_lal_1(xi_za_lal_1 <= 0) = [];
% xi_za_lal_1(xi_za_lal_1 <= 0) = [];
% f_za_hal_1(xi_za_hal_1 <= 0) = [];
% xi_za_hal_1(xi_za_hal_1 <= 0) = [];
% f_za_exdil(xi_za_exdil <= 0) = [];
% xi_za_exdil(xi_za_exdil <= 0) = [];
% f_za_excol(xi_za_excol <= 0) = [];
% xi_za_excol(xi_za_excol <= 0) = [];
% 
% scale33 = -0.04;
% 
% plot(scale33 * f_za_lal_1 + 0.7, xi_za_lal_1, 'Color', hex2rgb('#8D493A'),...
%     'LineWidth', 1.25)
% fill([scale33 * f_za_lal_1 + 0.7, 0.7 * ones(size(f_za_lal_1))],...
%      [xi_za_lal_1, fliplr(xi_za_lal_1)], hex2rgb('#DC8686'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale33 * f_za_hal_1 + 1.7, xi_za_hal_1, 'Color', hex2rgb('#537188'),...
%     'LineWidth', 1.25)
% fill([scale33 * f_za_hal_1 + 1.7, 1.7 * ones(size(f_za_hal_1))],...
%      [xi_za_hal_1, fliplr(xi_za_hal_1)], hex2rgb('#7EACB5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale33 * f_za_exdil + 2.7, xi_za_exdil, 'Color', hex2rgb('#8174A0'),...
%     'LineWidth', 1.25)
% fill([scale33 * f_za_exdil + 2.7, 2.7 * ones(size(f_za_exdil))],...
%      [xi_za_exdil, fliplr(xi_za_exdil)], hex2rgb('#A888B5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale33 * f_za_excol + 3.7, xi_za_excol, 'Color', hex2rgb('#659287'),...
%     'LineWidth', 1.25)
% fill([scale33 * f_za_excol + 3.7, 3.7 * ones(size(f_za_excol))],...
%      [xi_za_excol, fliplr(xi_za_excol)], hex2rgb('#B1C29E'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.02 0.02])
ylabel('$z_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
% xlim([0.3, 4.3])
ylim([0.1, 0.71])

%% optical sharpness distribution

nexttile

bp34 = boxplot([cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).sbar_opt);...
    cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).sbar_opt);...
    cat(1, Aggs_exdil.sbar_opt); cat(1, Aggs_excol.sbar_opt)],...
    condition31, 'Notch', 'on', 'Symbol', 'o', 'Widths', 0.25);

boxes34 = findobj(bp34, 'Tag', 'Box');
patch(get(boxes34(1), 'XData'), get(boxes34(1), 'YData'), hex2rgb('#DC8686'),...
    'EdgeColor', hex2rgb('#8D493A'), 'FaceAlpha', 0.5);
patch(get(boxes34(2), 'XData'), get(boxes34(2), 'YData'), hex2rgb('#7EACB5'),...
    'EdgeColor', hex2rgb('#537188'), 'FaceAlpha', 0.5);
patch(get(boxes34(3), 'XData'), get(boxes34(3), 'YData'), hex2rgb('#A888B5'),...
    'EdgeColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.5);
patch(get(boxes34(4), 'XData'), get(boxes34(4), 'YData'), hex2rgb('#B1C29E'),...
    'EdgeColor', hex2rgb('#659287'), 'FaceAlpha', 0.5);

medians34 = findobj(bp34, 'Tag', 'Median');
set(medians34(1), 'Color', hex2rgb('#632626'), 'LineWidth', 2);
set(medians34(2), 'Color', hex2rgb('#374259'), 'LineWidth', 2);
set(medians34(3), 'Color', hex2rgb('#574964'), 'LineWidth', 2);
set(medians34(4), 'Color', hex2rgb('#5F6F65'), 'LineWidth', 2);

outliers34 = findobj(bp34, 'Tag', 'Outliers');
outliers34(1).MarkerEdgeColor = hex2rgb('#DC8686');
outliers34(1).MarkerSize = 3;
outliers34(2).MarkerEdgeColor = hex2rgb('#7EACB5');
outliers34(2).MarkerSize = 3;
outliers34(3).MarkerEdgeColor = hex2rgb('#A888B5');
outliers34(3).MarkerSize = 3;
outliers34(4).MarkerEdgeColor = hex2rgb('#B1C29E');
outliers34(4).MarkerSize = 3;

upwhisker34 = findobj(gca,'type', 'line', 'tag', 'Upper Whisker');
set(upwhisker34, 'linestyle', '-');
lowwhisker34 = findobj(gca, 'type', 'line','tag', 'Lower Whisker');
set(lowwhisker34, 'linestyle', '-');

% hold on
% 
% % Compute kernel density estimate for each condition
% [f_sa_lal_1, xi_sa_lal_1] =...
%     ksdensity(cat(1, Aggs_lal_1(cat(1, Aggs_lal_1.n_subagg) > 0).sbar_opt));
% [f_sa_hal_1, xi_sa_hal_1] =...
%     ksdensity(cat(1, Aggs_hal_1(cat(1, Aggs_hal_1.n_subagg) > 0).sbar_opt));
% [f_sa_exdil, xi_sa_exdil] =...
%     ksdensity(cat(1, Aggs_exdil.sbar_opt));
% [f_sa_excol, xi_sa_excol] =...
%     ksdensity(cat(1, Aggs_excol.sbar_opt));
% 
% % to avoid issue with log scale in y axis
% f_sa_lal_1(xi_sa_lal_1 <= 0) = [];
% xi_sa_lal_1(xi_sa_lal_1 <= 0) = [];
% f_sa_hal_1(xi_sa_hal_1 <= 0) = [];
% xi_sa_hal_1(xi_sa_hal_1 <= 0) = [];
% f_sa_exdil(xi_sa_exdil <= 0) = [];
% xi_sa_exdil(xi_sa_exdil <= 0) = [];
% f_sa_excol(xi_sa_excol <= 0) = [];
% xi_sa_excol(xi_sa_excol <= 0) = [];
% 
% scale34 = -0.12;
% 
% plot(scale34 * f_sa_lal_1 + 0.7, xi_sa_lal_1, 'Color', hex2rgb('#8D493A'),...
%     'LineWidth', 1.25)
% fill([scale34 * f_sa_lal_1 + 0.7, 0.7 * ones(size(f_sa_lal_1))],...
%      [xi_sa_lal_1, fliplr(xi_sa_lal_1)], hex2rgb('#DC8686'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale34 * f_sa_hal_1 + 1.7, xi_sa_hal_1, 'Color', hex2rgb('#537188'),...
%     'LineWidth', 1.25)
% fill([scale34 * f_sa_hal_1 + 1.7, 1.7 * ones(size(f_sa_hal_1))],...
%      [xi_sa_hal_1, fliplr(xi_sa_hal_1)], hex2rgb('#7EACB5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale34 * f_sa_exdil + 2.7, xi_sa_exdil, 'Color', hex2rgb('#8174A0'),...
%     'LineWidth', 1.25)
% fill([scale34 * f_sa_exdil + 2.7, 2.7 * ones(size(f_sa_exdil))],...
%      [xi_sa_exdil, fliplr(xi_sa_exdil)], hex2rgb('#A888B5'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% plot(scale34 * f_sa_excol + 3.7, xi_sa_excol, 'Color', hex2rgb('#659287'),...
%     'LineWidth', 1.25)
% fill([scale34 * f_sa_excol + 3.7, 3.7 * ones(size(f_sa_excol))],...
%      [xi_sa_excol, fliplr(xi_sa_excol)], hex2rgb('#B1C29E'),...
%      'FaceAlpha', 0.5, 'EdgeColor', 'none');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.02 0.02])
ylabel('$s_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
% xlim([0.3, 4.3])
ylim([0.5, 2.75])

%% n_hyb freqeuncy comparison subplot %%

% initialize figure 3
f4 = figure;
f4.Position = [200, 150, 1200, 600];
set(f4, 'color', 'white');

n_subagg_0 = {cat(1,Aggs_lal_1.n_subagg), cat(1,Aggs_hal_1.n_subagg),...
    cat(1,Aggs_exdil.n_subagg), cat(1,Aggs_excol.n_subagg)};

tt4 = tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

xlbl4 = [strcat('Lo-Aglom', {' '}, '(n =', {' '}, num2str(n_aggs_tot(1)), ')'),...
    strcat('Mod-Colaps', {' '},  '(n =', {' '}, num2str(n_aggs_tot(2)), ')'),...
    strcat('Hi-Aglom', {' '},  '(n =', {' '}, num2str(length(n_subagg_0{3})), ')'),...
    strcat('Ex-Colaps', {' '},  '(n =', {' '}, num2str(length(n_subagg_0{4})), ')')];
condition4 = categorical(xlbl4, {xlbl4{1}, xlbl4{2}, xlbl4{3}, xlbl4{4}});

nexttile

n_subagg = n_subagg_0;

n_subagg{1}(n_subagg_0{1} < 1) = [];
n_subagg{2}(n_subagg_0{2} < 1) = [];

n_hyb = 100 * [nnz(n_subagg{1} == 1),...
    nnz(n_subagg{1} == 2),...
    nnz((n_subagg{1} >= 3) & (n_subagg{1} <= 5)),...
    nnz((n_subagg{1} >= 6) & (n_subagg{1} <= 10)),...
    nnz(n_subagg{1} > 10);...
    nnz(n_subagg{2} == 1),...
    nnz(n_subagg{2} == 2),...
    nnz((n_subagg{2} >= 3) & (n_subagg{2} <= 5)),...
    nnz((n_subagg{2} >= 6) & (n_subagg{2} <= 10)),...
    nnz(n_subagg{2} > 10);...
    nnz(n_subagg{3} == 1),...
    nnz(n_subagg{3} == 2),...
    nnz((n_subagg{3} >= 3) & (n_subagg{3} <= 5)),...
    nnz((n_subagg{3} >= 6) & (n_subagg{3} <= 10)),...
    nnz(n_subagg{3} > 10);...
    nnz(n_subagg{4} == 1),...
    nnz(n_subagg{4} == 2),...
    nnz((n_subagg{4} >= 3) & (n_subagg{4} <= 5)),...
    nnz((n_subagg{4} >= 6) & (n_subagg{4} <= 10)),...
    nnz(n_subagg{4} > 10)];
n_hyb(1,:) = n_hyb(1,:) / length(n_subagg{1});
n_hyb(2,:) = n_hyb(2,:) / length(n_subagg{2});
n_hyb(3,:) = n_hyb(3,:) / length(n_subagg{3});
n_hyb(4,:) = n_hyb(4,:) / length(n_subagg{4});

b41 = bar(condition4, n_hyb,'stacked');
clr41 = hex2rgb({'#C96868', '#FADFA1', '#7EACB5', '#8174A0', '#574964'});

for i = 1 : size(n_hyb,2)
    b41(i).BarWidth = 0.4;
    b41(i).FaceColor = clr41(i,:);
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02])
ylabel('Frequency [$\%$]', 'interpreter', 'latex', 'FontSize', 22)
yticks([0 20 40 60 80 100])
ylim([0 100])

lgd41 = legend({'$n_\mathrm{hyb} = 1$', '$n_\mathrm{hyb} = 2$',...
    '$3 \le n_\mathrm{hyb} \le 5$', '$6 \le n_\mathrm{hyb} \le 10$',...
    '$n_\mathrm{hyb} > 10$'}, 'interpreter', 'latex', 'FontSize', 18,...
    'location', 'northoutside', 'orientation', 'horizontal',...
    'NumColumns', 2);
lgd41.ItemTokenSize = [15, 15];

% manually positioned downward arrows on n_hyb = 1 in left panel
% adjust these to visually align with bar groups

xBarPos = [1.13, 2.04, 3, 4];  % tweak as needed for alignment

% axis and position
ax4_1 = gca;
axPos = ax4_1.Position;
ylim4 = ax4_1.YLim;

for groupIdx = 1:2
    xBar = xBarPos(groupIdx);
    yBarTop = b41(1).YData(groupIdx);  % top of red bar

    % normalized figure coordinates
    xFig = axPos(1) + (xBar - 0.5)/4 * axPos(3);  % assuming xlim = [0.5, 4.5]
    yFig = axPos(2) + (yBarTop - ylim4(1)) / diff(ylim4) * axPos(4);

    dy = 0.05;  % arrow length in normalized units

    annotation('textarrow', [xFig xFig], [yFig yFig - dy], ...
        'String', '', ...
        'HeadStyle', 'plain', ...        
        'LineWidth', 2.0, ...            
        'Color', [0 0 0]);               
end

nexttile

n_colaps = {cat(1,Aggs_lal_1.n_colaps), cat(1,Aggs_hal_1.n_colaps),...
    cat(1,Aggs_exdil.n_colaps), cat(1,Aggs_excol.n_colaps)};

n_colaps{1}(n_subagg_0{1} < 1) = [];
r0_colaps{1} = double(n_colaps{1}) ./ double(n_subagg{1});
n_colaps{2}(n_subagg_0{2} < 1) = [];
r0_colaps{2} = double(n_colaps{2}) ./ double(n_subagg{2});
r0_colaps{3} = double(n_colaps{3}) ./ double(n_subagg{3});
r0_colaps{4} = double(n_colaps{4}) ./ double(n_subagg{4});

r_colaps = 100 * [nnz(r0_colaps{1} == 0),...
    nnz((r0_colaps{1} > 0) & (r0_colaps{1} < 0.33)),...
    nnz((r0_colaps{1} >= 0.33) & (r0_colaps{1} <= 0.67)),...
    nnz((r0_colaps{1} > 0.67) & (r0_colaps{1} < 1)),...
    nnz(r0_colaps{1} == 1);...
    nnz(r0_colaps{2} == 0),...
    nnz((r0_colaps{2} > 0) & (r0_colaps{2} < 0.33)),...
    nnz((r0_colaps{2} >= 0.33) & (r0_colaps{2} <= 0.67)),...
    nnz((r0_colaps{2} > 0.67) & (r0_colaps{2} < 1)),...
    nnz(r0_colaps{2} == 1);...
    nnz(r0_colaps{3} == 0),...
    nnz((r0_colaps{3} > 0) & (r0_colaps{3} < 0.33)),...
    nnz((r0_colaps{3} >= 0.33) & (r0_colaps{3} <= 0.67)),...
    nnz((r0_colaps{3} > 0.67) & (r0_colaps{3} < 1)),...
    nnz(r0_colaps{3} == 1);...
    nnz(r0_colaps{4} == 0),...
    nnz((r0_colaps{4} > 0) & (r0_colaps{4} < 0.33)),...
    nnz((r0_colaps{4} >= 0.33) & (r0_colaps{4} <= 0.67)),...
    nnz((r0_colaps{4} > 0.67) & (r0_colaps{4} < 1)),...
    nnz(r0_colaps{4} == 1)];
r_colaps(1,:) = r_colaps(1,:) / length(n_subagg{1});
r_colaps(2,:) = r_colaps(2,:) / length(n_subagg{2});
r_colaps(3,:) = r_colaps(3,:) / length(n_subagg{3});
r_colaps(4,:) = r_colaps(4,:) / length(n_subagg{4});

condition4 = categorical(xlbl4, {xlbl4{1}, xlbl4{2}, xlbl4{3},...
    xlbl4{4}});

b42 = bar(condition4, r_colaps,'stacked');
clr42 = hex2rgb({'#DEAA79', '#FFE6A9', '#B1C29E', '#659287', '#3A4D39'});

for i = 1 : size(r_colaps,2)
    b42(i).BarWidth = 0.4;
    b42(i).FaceColor = clr42(i,:);
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02])
yticks([0 20 40 60 80 100])
ylim([0 100])

lgd42 = legend({'$r_\mathrm{clp} = 0$', '$0 < r_\mathrm{clp} < 0.33$',...
    '$0.33 \le r_\mathrm{clp} \le 0.67$', '$0.67 < r_\mathrm{clp} < 1$',...
    '$r_\mathrm{clp} = 1$'}, 'interpreter', 'latex', 'FontSize', 18,...
    'location', 'northoutside', 'orientation', 'horizontal', 'NumColumns', 2);
lgd42.ItemTokenSize = [15, 15];

%% dpp vs. da as a function of n_hyb, not experimental case %%

% initialize figure 5
f5 = figure;
f5.Position = [250, 200, 550, 600];
set(f5, 'color', 'white')

plt5 = cell(5,1); % initialize dpp vs da plots per n_hyb

 % compile number of subaggregates
n_subagg_tot = cat(1, n_subagg_0{1}(id_agg_lal_1), n_subagg_0{2}(id_agg_hal_1));

% compile mean primary particle sizes
dbarpp_tot = [dbarpp_manu_lal_1; dbarpp_manu_hal_1];

% compile projected-area sizes
da_tot = [cat(1, Aggs_lal_1(id_agg_lal_1).da);...
    cat(1, Aggs_hal_1(id_agg_hal_1).da)];

plt5{1} = plot(da0, dpp0, 'Color', hex2rgb('#DEAA79'),... % [0.4940 0.1840 0.5560]
    'LineStyle', '-.', 'LineWidth', 2.5);
hold on

plt5{2} = scatter(da_tot(n_subagg_tot == 1), dbarpp_tot(n_subagg_tot == 1),...
    25, clr41(1,:), '^', 'LineWidth', 1.5);

plt5{3} = scatter(da_tot(n_subagg_tot == 2), dbarpp_tot(n_subagg_tot == 2),...
    35, clr41(2,:), 's', 'LineWidth', 1.5);

plt5{4} = scatter(da_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)),...
    dbarpp_tot((n_subagg_tot >= 3) & (n_subagg_tot <= 5)), 35,...
    clr41(3,:), 'h', 'LineWidth', 1.5);

plt5{5} = scatter(da_tot(n_subagg_tot > 5), dbarpp_tot(n_subagg_tot > 5),...
    25, clr41(4,:), 'o', 'LineWidth', 1.5);

% plot configs
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * xmin1, 1.2 * xmax1])
ylim([0.95 * ymin1, 1.05 * ymax1])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\overline{d}_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

% number of aggregates to be manually sized
n_aggs_manu = [size(dpp_manu_lal_1, 1), size(dpp_manu_hal_1, 1)];

legend(cat(2, plt5{:}),...
    {'Olfert and Rogak (2019)', '$n_\mathrm{hyb} = 1$',...
    '$n_\mathrm{hyb} = 2$', '$3 \le n_\mathrm{hyb} \le 5$',...
    '$n_\mathrm{hyb} > 5$'}, 'interpreter', 'latex', 'FontSize', 12,...
    'location', 'northoutside', 'orientation', 'horizontal',...
    'NumColumns', 2)

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
