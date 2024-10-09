clc
clear
close all
warning('off')

%% initialize dpp vs. da figure %%

f1 = figure;
f1.Position = [100, 100, 750, 500];
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

fname_agg_lal_1 = '19AUG24_LAL_End_Slider.mat';
fdir_agg_lal_1 = 'D:\HN\AUG24Onward\TEM\SimMag-07OCT-24-0830pm\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ATEMS_Area';
fname_pp_lal_1 = 'PFA_ET+NIT_LAL_19AUG24_End';
fdir_pp_lal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ImageJ_Primaries\CSV';
id_agg_lal_1 = [39, 42, 44, 56, 57, 59, 58, 62, 63, 64];

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

n_agg_lal_1 = length(id_agg_lal_1);
dpp_manu_lal_1 = cell(n_agg_lal_1, 1);
dbarpp_manu_lal_1 = zeros(n_agg_lal_1, 1);
sigmapp_manu_lal_1 = zeros(n_agg_lal_1, 1);
npp_manu_lal_1 = zeros(n_agg_lal_1, 1);

fadd_pp_lal_1 = cell(n_agg_lal_1, 1);

figure(f1)
for i = 1 : n_agg_lal_1
    
    fadd_pp_lal_1{i} = char(strcat(fdir_pp_lal_1, {'\'}, fname_pp_lal_1,...
        '_', id_agg_lal_1(i), '.csv'));

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
fname_agg_hal_1 = '28AUG24_HAL_End_Slider.mat';
fdir_agg_hal_1 = 'D:\HN\AUG24Onward\TEM\SimMag-07OCT-24-0830pm\26SEP24_PFA_ET+NIT_HAL_28AUG24_End\ATEMS_Area';
fname_pp_hal_1 = 'PFA_ET+NIT_28AUG24_HAL_End';
fdir_pp_hal_1 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\26SEP24_PFA_ET+NIT_HAL_28AUG24_End\ImageJ_Primaries\CSV';
id_agg_hal_1 = 1:10;

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

n_agg_hal_1 = length(id_agg_hal_1);

dpp_manu_hal_1 = cell(n_agg_hal_1, 1);
dbarpp_manu_hal_1 = zeros(n_agg_hal_1, 1);
sigmapp_manu_hal_1 = zeros(n_agg_hal_1, 1);
npp_manu_hal_1 = zeros(n_agg_hal_1, 1);

fadd_pp_hal_1 = cell(n_agg_hal_1, 1);

for i = 1 : n_agg_hal_1

    fadd_pp_hal_1{i} = char(strcat(fdir_pp_hal_1, {'\'}, fname_pp_hal_1,...
        '_', id_agg_hal_1(i), '.csv'));

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

% %% "low" agglomeration level data - Repeat 2 %%
% 
% % load files for the aggregate segmentation data
% fname_agg_lal_2 = '20AUG-LAL-Start_Slider.mat';
% fdir_agg_lal_2 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ATEMS_Area_New';
% id_agg_lal_2 = [4, 5, 7, 11, 14, 15, 16, 27];
% fadd_agg_lal_2 = cell2mat(strcat(fdir_agg_lal_2, {'\'}, fname_agg_lal_2));
% load(fadd_agg_lal_2);
% 
% % Rename Aggs variable
% 
% vars = who; % call variable names
% 
% for i = 1 : length(vars)
% 
%     varname = vars{i};  % Get the name of the variable as a string
% 
%     if strcmp(varname, 'Aggs') % check if the variable contains Aggs
%         newVarName = ['Aggs' '_lal_2']; % Modify the variable name by adding a suffix
% 
%         % Assign the value of the original variable to the new variable...
%         %   ...(i.e. dynamically create a new variable)
%         eval([newVarName ' = ' varname ';']);
%         clear(varname); % Delete the old variable
%     end
% 
% end
% 
% clear vars varname newVarName
% 
% % information on primary particle manual sizing data
% fname_pp_lal_2 = {'_20AUG24_LAL_4', '_20AUG24_LAL_5',...
%     '_20AUG24_LAL_7', '_20AUG24_LAL_11',...
%     '_20AUG24_LAL_14', '_20AUG24_LAL_15',...
%     '_20AUG24_LAL_16', '_20AUG24_LAL_27'};
% fdir_pp_lal_2 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ImageJ_Primaries_R2\20AUG24_LAL_Start\CSV';
% 
% n_agg_lal_2 = length(id_agg_lal_2);
% 
% dpp_manu_lal_2 = cell(n_agg_lal_2, 1);
% dbarpp_manu_lal_2 = zeros(n_agg_lal_2, 1);
% sigmapp_manu_lal_2 = zeros(n_agg_lal_2, 1);
% npp_manu_lal_2 = zeros(n_agg_lal_2, 1);
% 
% fadd_pp_lal_2 = cell(n_agg_lal_2, 1);
% 
% for i = 1 : n_agg_lal_2
% 
%     fadd_pp_lal_2{i} = char(strcat(fdir_pp_lal_2, {'\'}, fname_pp_lal_2{i}, '.csv'));
% 
%     if isfile(fadd_pp_lal_2{i})
% 
%         % opts_pp = detectImportOptions(fname_pp{i});
%         % opts_pp = setvartype(opts_pp, 'char');
%         pp_manu_lal_2 = readtable(fadd_pp_lal_2{i});
% 
%         dpp_manu_lal_2{i} = sqrt(4 * pp_manu_lal_2.Area(2:end) / pi);
% 
%         dbarpp_manu_lal_2(i) = geomean(dpp_manu_lal_2{i});
%         sigmapp_manu_lal_2(i) = morph.geostd(dpp_manu_lal_2{i});
%         npp_manu_lal_2(i) = length(dpp_manu_lal_2{i});
% 
%         plt_lal_2 = scatter(Aggs_lal_2(id_agg_lal_2(i)).da, dbarpp_manu_lal_2(i), 15,...
%             hex2rgb('#006989'), 'v');
% 
%     end
% 
% end
% 
% %% "high" agglomeration level data - Repeat 2 %%
% 
% % load files for the aggregate segmentation data
% fname_agg_hal_2 = '28AUG-HAL-Start_Slider.mat';
% fdir_agg_hal_2 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ATEMS_Area_New';
% id_agg_hal_2 = [9, 10, 13, 20, 21, 22, 31, 43, 50, 82];
% fadd_agg_hal_2 = cell2mat(strcat(fdir_agg_hal_2, {'\'}, fname_agg_hal_2));
% load(fadd_agg_hal_2);
% 
% % Rename Aggs variable
% 
% vars = who; % call variable names
% 
% for i = 1 : length(vars)
% 
%     varname = vars{i};  % Get the name of the variable as a string
% 
%     if strcmp(varname, 'Aggs') % check if the variable contains Aggs
%         newVarName = ['Aggs' '_hal_2']; % Modify the variable name by adding a suffix
% 
%         % Assign the value of the original variable to the new variable...
%         %   ...(i.e. dynamically create a new variable)
%         eval([newVarName ' = ' varname ';']);
%         clear(varname); % Delete the old variable
%     end
% 
% end
% 
% clear vars varname newVarName
% 
% % information on primary particle manual sizing data
% fname_pp_hal_2 = {'_28AUG24_HAL_9', '_28AUG24_HAL_10',...
%     '_28AUG24_HAL_13', '_28AUG24_HAL_20',...
%     '_28AUG24_HAL_21', '_28AUG24_HAL_22',...
%     '_28AUG24_HAL_31', '_28AUG24_HAL_43',...
%     '_28AUG24_HAL_50', '_28AUG24_HAL_82'};
% fdir_pp_hal_2 = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\ImageJ_Primaries_R2\28AUG24_HAL_Start\CSV';
% 
% n_agg_hal_2 = length(id_agg_hal_2);
% 
% dpp_manu_hal_2 = cell(n_agg_hal_2, 1);
% dbarpp_manu_hal_2 = zeros(n_agg_hal_2, 1);
% sigmapp_manu_hal_2 = zeros(n_agg_hal_2, 1);
% npp_manu_hal_2 = zeros(n_agg_hal_2, 1);
% 
% fadd_pp_hal_2 = cell(n_agg_hal_2, 1);
% 
% for i = 1 : n_agg_hal_2
% 
%     fadd_pp_hal_2{i} = char(strcat(fdir_pp_hal_2, {'\'}, fname_pp_hal_2{i}, '.csv'));
% 
%     if isfile(fadd_pp_hal_2{i})
% 
%         % opts_pp = detectImportOptions(fname_pp{i});
%         % opts_pp = setvartype(opts_pp, 'char');
%         pp_manu_hal_2 = readtable(fadd_pp_hal_2{i});
% 
%         dpp_manu_hal_2{i} = sqrt(4 * pp_manu_hal_2.Area(2:end) / pi);
% 
%         dbarpp_manu_hal_2(i) = geomean(dpp_manu_hal_2{i});
%         sigmapp_manu_hal_2(i) = morph.geostd(dpp_manu_hal_2{i});
%         npp_manu_hal_2(i) = length(dpp_manu_hal_2{i});
% 
%         plt_hal_2 = scatter(Aggs_hal_2(id_agg_hal_2(i)).da, dbarpp_manu_hal_2(i), 15,...
%             hex2rgb('#C96868'), 'd');
% 
%     end
% 
% end
% 
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([10,1500])
ylim([5,50])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

% legend(cat(2, plt_0, plt_lal_1, plt_lal_2, plt_hal_1, plt_hal_2),...
%     cat(2, {'Olfert and Rogak (2019)'},...
%     {'Low agglom. - Rep. 1'}, {'Low agglom. - Rep. 2'},...
%     {'High agglom. - Rep. 1'}, {'High agglom. - Rep. 2'}),...
%     'interpreter', 'latex', 'FontSize', 14, 'location', 'northeastoutside')
legend(cat(2, plt_0, plt_lal_1, plt_hal_1),...
    cat(2, {'Olfert and Rogak (2019)'},...
    {'Low agglom.'}, {'High agglom.'}),...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'northeastoutside')


%% da comparison figure %%

f2 = figure;
f2.Position = [100, 100, 500, 500];
set(f2, 'color', 'white');

% condition = categorical([repmat({'Low agglom. - Rep. 1'}, size(Aggs_lal_1, 2), 1);...
%     repmat({'Low agglom. - Rep. 2'}, size(Aggs_lal_2, 2), 1);...
%     repmat({'High agglom. - Rep. 1'}, size(Aggs_hal_1, 2), 1);...
%     repmat({'High agglom. - Rep. 2'}, size(Aggs_hal_2, 2), 1)]);
condition2 = categorical([repmat({'Low agglom.'}, size(Aggs_lal_1, 2), 1);...
    repmat({'High agglom.'}, size(Aggs_hal_1, 2), 1)]);

% boxplot([cat(1, Aggs_lal_1.da); cat(1, Aggs_lal_2.da);...
%     cat(1, Aggs_hal_1.da); cat(1, Aggs_hal_2.da)], condition)
boxplot([cat(1, Aggs_lal_1.da); cat(1, Aggs_hal_1.da)], condition2)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
% ylim([10,60])
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% dpp_agg comparison figure %%

f3 = figure;
f3.Position = [100, 100, 500, 500];
set(f3, 'color', 'white');

condition3 = categorical([repmat({'Low agglom.'}, size(dpp_manu_lal_1, 1), 1);...
    repmat({'High agglom.'}, size(dpp_manu_hal_1, 1), 1)]);

boxplot([dbarpp_manu_lal_1; dbarpp_manu_hal_1], condition3)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
% ylim([10,60])
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp_{agg}}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)


%% dpp_ens comparison figure %%

f4 = figure;
f4.Position = [100, 100, 500, 500];
set(f4, 'color', 'white');

dpp_ens_lal_1 = cat(1, cell2mat(dpp_manu_lal_1));
dpp_ens_hal_1 = cat(1, cell2mat(dpp_manu_hal_1));

condition4 = categorical([repmat({'Low agglom.'}, size(dpp_ens_lal_1, 1), 1);...
    repmat({'High agglom.'}, size(dpp_ens_hal_1, 1), 1)]);

boxplot([dpp_ens_lal_1; dpp_ens_hal_1], condition4)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10,...
    'TickLength', [0.02 0.02], 'YScale', 'log')
% ylim([10,60])
xlabel('Condition', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp_{ens}}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
