
% SCATTERGRAPH_UI   Draws a scatter plot comparing 2 characteristics of an 
%                   Aggs struct. 
% Author:           Darwin Zhu, 2021-07-01

%=========================================================================%

function [] = scattergraph_ui()
clear all;
close all;
clc;

hold on;
% Run load_all to compile all case files into struct h.

f_finished = 0;
legendlabels = {};
indx_x = 1;
indx_y = 1;
f_firstrun = 0;

while f_finished == 0
    uiwait(msgbox('Please navigate to the .mat file containing the Aggs_pp struct',...
            'Loading aggregate struct...','help'));
    [baseFileName, folder] = uigetfile('*.mat');
    fullFileName = fullfile(folder, baseFileName);

    if exist(fullFileName, 'file')
        % Normal situation - they picked an existing file.
        file_in = load(fullFileName);
        aggname = fieldnames(file_in);
        Aggs = file_in.(aggname{1});
        % Now do something with storedStructure, like extract fields into new variables or whatever you want.
    else
        % Error: Would only get here if they typed in a name of a non-existant file
        % instead of picking one from the folder.
        warningMessage = sprintf('Warning: mat file does not exist:\n%s', fullFileName);
        uiwait(errordlg(warningMessage));
        return;
    end

    % Get all the fields.
    fields = fieldnames(Aggs);

    if (f_firstrun == 0) 
        % Ask for which 2 properties to compare.
        [indx_x,tf_x] = listdlg('PromptString',{'Select a property to plot on the x-axis.',...
            'Only one file can be selected at a time.',''},...
            'SelectionMode','single','ListString',fields);

        [indx_y,tf_y] = listdlg('PromptString',{'Select a property to plot on the y-axis.',...
            'Only one file can be selected at a time.',''},...
            'SelectionMode','single','ListString',fields);
        f_firstrun = 1;
    end

    if (tf_x == 0 || tf_y == 0) 
        uiwait('Error: property left blank');
        return;
    end

    property_x = [Aggs.(fields{indx_x})];
    property_y = [Aggs.(fields{indx_y})];

    % Create an errorbar plot. 
     scatter(property_x, property_y, [],[rand, rand, rand], 'o');
     xlabel(fields{indx_x});
     ylabel(fields{indx_y});
     
     
     legendlabels{end+1} = baseFileName(1:end-4);
     
     choice = questdlg('Do you want to plot another aggregate structure ?',...
        'Continue?', 'Yes', 'No', 'Yes');
        if strcmp(choice,'Yes')
        	f_finished = 0;
        else
        	f_finished = 1;
        end
     
end

%Plot legend without LaTeX interpreter messing with underscores.
legend(legendlabels, 'Interpreter', 'none');
hold off;

end