
% HDF_PLOT          Draws a scatter plot of two characteristics of a given
%                   HDF5 file. User selects the properties to plot, whether
%                   to set x and y axes to log scale, and number of files
%                   to plot.
%
% INPUTS: NONE
% OUTPUTS: NONE
%
% Author:           Darwin Zhu, 2021-08-10

%=========================================================================%

function [] = hdf_plot()
    clear;
    close all;
    clc;
    
    hold on;
    
    % Initialize fields
    legendlabels = {};
    fields = {};
    f_xlog = 0;
    f_ylog = 0;
    f_finished = 0;
    f_firstrun = 0;
    h5names = {};
    monitoring_select = 0;
    
    % Random color setup
    r = 0.8;
    g = 0;
    b = 0;
    
    
    % Preserve the original location
    start_path = pwd;
    
    % Loop until user exits out of loop
    while f_finished == 0
        
        % Select file to plot.
        uiwait(msgbox('Please navigate to the .h5 file to plot.',...
                    'Loading files...','help'));
        [baseFileName, folder] = uigetfile('*.h5')
        cd(folder);
        h5metadata = h5info(baseFileName);

        
        % Only select the x-attribute once; we assume future attributes are
        % plotted against the same axis.
        if (f_firstrun == 0)
           
            choice = questdlg('Plot MonitoringData or ParticleData?',...
            'Directory selection', 'MonitoringData', 'ParticleData', 'File');
            if strcmp(choice,'MonitoringData')
                h5datasets = h5metadata.Groups.Groups(1).Datasets;
                monitoring_select = 1;
            else
                h5datasets = h5metadata.Groups.Groups(2).Datasets;
                monitoring_select = 0;
            end

            % Iterate over h5datasets.Name, adding items to list of names
            dim = size(h5datasets);
            for i = 1:dim(1)
                fields{end+1} = h5datasets(i).Name;
            end
        
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
        end
        
         
        
        % Only ask once at the beginning if y axis should be in log scale.
        if (f_firstrun == 0)
            choice = questdlg('Set y axis to log scale?',...
            'Continue?', 'Yes', 'No', 'Yes');
            if strcmp(choice,'Yes')
                f_ylog = 1;
            else
                f_ylog = 0;
            end
        end
        
        % Exit function if either property is left blank.
        if (tf_x == 0 || tf_y == 0) 
            uiwait('Error: property left blank');
            return;
        end
        
        % Read data from correct data path depending on prior choice.
        if(monitoring_select == 0) 
            property_x = h5read(baseFileName, append('/NEO/ParticleData/', fields{indx_x}));
            
            % Check if they selected a FluorPeak option as the x axis.
            if strcmp(fields{indx_x}, 'Xe1_FluorPeak') || strcmp(fields{indx_x}, 'Xe2_FluorPeak')
                property_x = property_x(2, :);
                property_x1 = property_x(1,:);
            end
            
            property_y = h5read(baseFileName, append('/NEO/ParticleData/', fields{indx_y}));
            % Check if they selected a FluorPeak option as the y axis.
            if strcmp(fields{indx_y}, 'Xe1_FluorPeak') || strcmp(fields{indx_y}, 'Xe2_FluorPeak')
                property_y = property_y(2, :);
                property_y1 = property_y(1,:);
            end
        else
            property_x = h5read(baseFileName, append('/NEO/MonitoringData/', fields{indx_x}));
            
            property_y = h5read(baseFileName, append('/NEO/MonitoringData/', fields{indx_y}));
            
        end
        
        % Create an errorbar plot. 
        scatter(property_x, property_y, [],[r, g, b], 'x');
        
        % Generate a new color
        [r,g,b] = randomColor(r,g,b);
        
        % Plot the extra Xe1_Fluorpeak column if it was selected.
        if strcmp(fields{indx_x}, 'Xe1_FluorPeak')
            scatter(property_x1, property_y, [],[rand, rand, rand], 'x');
        elseif strcmp(fields{indx_y}, 'Xe1_FluorPeak')
            scatter(property_x, property_y1, [],[rand, rand, rand], 'x');
        end
        xlabel({fields{indx_x}}, 'Interpreter', 'none');
        ylabel({fields{indx_y}}, 'Interpreter', 'none');
        
        % Set log scales based on user response
        if (f_xlog == 1) 
            set(gca, 'XScale', 'log')
        end
        if (f_ylog == 1)
            set(gca, 'YScale', 'log')
        end
        
        % Add name to legends list.
        if strcmp(fields{indx_x}, 'Xe1_FluorPeak') || strcmp(fields{indx_y}, 'Xe1_FluorPeak')
            legendlabels{end+1} = append(fields{indx_y}, ', Xe1_Fluor Ch2');
            legendlabels{end+1} = append(fields{indx_y}, ', Xe1_Fluor Ch1');
        else
            legendlabels{end+1} = fields{indx_y};
        end
        
        % Check if user wants to add an additional property to plot.
        choice = questdlg('Do you want to plot another file?',...
        'Continue?', 'Yes', 'No', 'Yes');
        if strcmp(choice,'Yes')
        	f_finished = 0;
        else
        	f_finished = 1;
        end
        
        % Clear field names.
        %fields = {};
        
        % Indicate that we have run through loop once.
        f_firstrun = 1;
    end
    
    %Plot legend without LaTeX interpreter messing with underscores.
    legend(legendlabels, 'Interpreter', 'none');
    hold off;
    
    % Reset to original file location
    cd(pwd);
end

% Function that converts rgb to hsl, changes hue, then switches back to
% rgb for output.
function [r2,g2,b2] = randomColor(r1,g1,b1)

    max_rgb = max([r1,g1,b1]);
    min_rgb =min([r1,g1,b1]);
    lum = (min_rgb + max_rgb)/2;
    if (lum <= 0.5)
        sat = (max_rgb - min_rgb)/(max_rgb + min_rgb);
    else
        sat = (max_rgb - min_rgb)/(2.0 - max_rgb - min_rgb);
    end
    
    if (max_rgb == r1)
        hue = (g1-b1)/(max_rgb-min_rgb);
    elseif (max_rgb == g1)
        hue = 2.0+(b1-r1)/(max_rgb-min_rgb);
    else
        hue = 4.0+(r1-g1)/(max_rgb-min_rgb);
    end
    
    % Shift hue value by a set value
    hue = hue * 60 + 137.508;
    if hue >= 360
        hue = hue - 360;
    end
    hue = hue/360;
    
    % Convert back from hsl to rgb
    if sat == 0
        r2 = lum;
        g2 = lum;
        b2 = lum;
    else
        if lum < 0.5
            temp_1 = lum*(1.0+sat);
        else
            temp_1 = lum + sat - lum*sat;
        end
        
        temp_2 = 2*lum-temp_1;
        
        temp_R = hue + 0.333;
        temp_G = hue;
        temp_B = hue - 0.333;
        r2 = threeTest(normalizeOne(temp_R), temp_1, temp_2);
        g2 = threeTest(normalizeOne(temp_G),temp_1, temp_2);
        b2 = threeTest(normalizeOne(temp_B),temp_1, temp_2);
    end
end

% Sets a value to be within the range 0-1.
function j = normalizeOne(i)
    while i < 0
        i = i + 1;
    end
    while i > 1
        i = i -1;
    end
    j = i;
end

function temp_out = threeTest(temp_x, temp_1, temp_2)
    if temp_x*6 < 1
        temp_out = temp_2 + (temp_1 - temp_2)*6*temp_x;
    elseif temp_x*2 < 1
        temp_out = temp_1;
    elseif temp_x*3 < 2
        temp_out = temp_2 + (temp_1 - temp_2)*(0.666-temp_x)*6;
    else
        temp_out = temp_2;
    end
end

