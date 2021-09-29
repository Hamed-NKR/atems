
% SCATTERGRAPH_UI   Draws a scatter plot comparing 2 characteristics of an 
%                   Aggs struct. Requires user to input which
%                   characteristics to plot, whether to plot multiple
%                   structs, and whether to select from file or folder. 
%
% INPUTS: NONE
% OUTPUTS: NONE
%
% Author:           Darwin Zhu, 2021-07-01

%=========================================================================%

function [] = scattergraph_ui()
clear all;
close all;
clc;

hold on;

% Initialize fields
f_finished = 0;
legendlabels = {};
indx_x = 1;
indx_y = 1;
f_firstrun = 0;
f_xlog = 0;
f_ylog = 0;
f_folder = 0;

% Loop while user still has aggregate files to plot.
while f_finished == 0
    
    % Choose either file selection or folder selection mode. 
    choice = questdlg('Select file or folder?',...
    'Directory selection', 'File', 'Folder', 'File');
    if strcmp(choice,'Folder')
        % f_folder indicates that folder option was chosen.
        f_folder = 1;
        
        % Open folder navigation UI in load_all and retrieve folder
        % contents.
        [baseFileName, h] = tools.load_all();
        
        % Recompile all rows in h into one Aggs struct.
        a1 = size(h);
        Aggs = struct([]);
        for i = 1:a1(2)
            Aggs = [Aggs, h(i).value];
        end
    else
        % f_folder indicates that folder option was not chosen.
        f_folder = 0;
        
        % Retrieve the file using uigetfile. Get the file name.
        uiwait(msgbox('Please navigate to the .mat file containing the Aggs_pp struct',...
                'Loading aggregate struct...','help'));
        [baseFileName, folder] = uigetfile('*.mat');
        fullFileName = fullfile(folder, baseFileName);

        % Check if file was correctly selected.
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
    end

    % Get all the fields.
    fields = fieldnames(Aggs);

    % Runs once at the beginning to select 2 properties to plot,, and assumes that subsequent aggs
    % structures will plot the same 2 properties.
    if (f_firstrun == 0) 
        % Ask for x-axis property from user.
        [indx_x,tf_x] = listdlg('PromptString',{'Select a property to plot on the x-axis.',...
            'Only one file can be selected at a time.',''},...
            'SelectionMode','single','ListString',fields);

        % Ask whether or not to set x axis to log scale.
        choice = questdlg('Set x axis to log scale?',...
        'Continue?', 'Yes', 'No', 'Yes');
        if strcmp(choice,'Yes')
        	f_xlog = 1;
        else
        	f_xlog = 0;
        end
        
        % Ask for y-axis property from user.
        [indx_y,tf_y] = listdlg('PromptString',{'Select a property to plot on the y-axis.',...
            'Only one file can be selected at a time.',''},...
            'SelectionMode','single','ListString',fields);
        
        % Ask whether or not to set y axis to log scale.
        choice = questdlg('Set y axis to log scale?',...
        'Continue?', 'Yes', 'No', 'Yes');
        if strcmp(choice,'Yes')
        	f_ylog = 1;
        else
        	f_ylog = 0;
        end
        
        % Indicate that loop has run once.
        f_firstrun = 1;
    end

    % 
    if (tf_x == 0 || tf_y == 0) 
        uiwait('Error: property left blank');
        return;
    end

    % Indicate the 2 properties' field names.
    property_x = [Aggs.(fields{indx_x})];
    property_y = [Aggs.(fields{indx_y})];

    % Create a scatterplot. 
     scatter(property_x, property_y, [],[rand, rand, rand], 'o');
     
     
    
     % Set log scales based on user response. Label if set to log scale.
     if (f_xlog == 1) 
         set(gca, 'XScale', 'log')
         xlabel(append('log(', fields{indx_x}, ')'));
     else
         xlabel(fields{indx_x});
     end
     
     if (f_ylog == 1)
         set(gca, 'YScale', 'log')
         ylabel(append('log(', fields{indx_y}, ')'));
     else
         ylabel(fields{indx_y});
     end
     
         
     
     if f_folder
         legendlabels{end+1} = baseFileName;
     else
         legendlabels{end+1} = baseFileName(1:end-4);
     end
     
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