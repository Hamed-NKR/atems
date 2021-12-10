% SCATTERGRAPH_DP_DA   Draws a scatter plot comparing manual and pcm dp to
%                      manual and kmeans da. User must navigate to the
%                      kmeans and manual sizing mat files. 
%
% INPUTS: NONE
% OUTPUTS: NONE
%
% Author:           Darwin Zhu, 2021-08-21

%=========================================================================%

function [Aggs_manual, Aggs_kmeans] = scattergraph_dp_da()
    clear all;
    close all;
    
    hold on;

    % Load in manual file using ui
    uiwait(msgbox('Please navigate to the folder containing the manual structs.',...
            'Loading aggregate struct...','help'));
        
    % Open folder navigation UI in load_all and retrieve folder
    % contents.
    [manual_file, h] = tools.load_all();

    % Recompile all rows in h into one Aggs struct.
    a1 = size(h);
    Aggs_manual = struct([]);
    for i = 1:a1(2)
        Aggs_manual = [Aggs_manual, h(i).value];
    end
    
    % Load kmeans file using ui
    uiwait(msgbox('Please navigate to the .mat file containing the kmeans struct',...
            'Loading aggregate struct...','help'));
    [kmeans_file, kmeans_folder] = uigetfile('*.mat');
    kmeans_full = fullfile(kmeans_folder, kmeans_file);
    
    % Check if file was correctly selected.
    if exist(kmeans_full, 'file')
        % Normal situation - they picked an existing file.
        kmeans_in = load(kmeans_full);
        kmeans_aggname = fieldnames(kmeans_in);
        Aggs_kmeans = kmeans_in.(kmeans_aggname{1});
        % Now do something with storedStructure, like extract fields into new variables or whatever you want.
    else
        % Error: Would only get here if they typed in a name of a non-existant file
        % instead of picking one from the folder.
        warningMessage = sprintf('Warning: mat file does not exist:\n%s', kmeans_full);
        uiwait(errordlg(warningMessage));
        return;
    end
    
    % Load carboseg file using ui
    uiwait(msgbox('Please navigate to the .mat file containing the carboseg struct',...
            'Loading aggregate struct...','help'));
    [cnn_file, cnn_folder] = uigetfile('*.mat');
    cnn_full = fullfile(cnn_folder, cnn_file);
    
    % Check if file was correctly selected.
    if exist(cnn_full, 'file')
        % Normal situation - they picked an existing file.
        cnn_in = load(cnn_full);
        cnn_aggname = fieldnames(cnn_in);
        Aggs_cnn = cnn_in.(cnn_aggname{1});
        % Now do something with storedStructure, like extract fields into new variables or whatever you want.
    else
        % Error: Would only get here if they typed in a name of a non-existant file
        % instead of picking one from the folder.
        warningMessage = sprintf('Warning: mat file does not exist:\n%s', cnn_full);
        uiwait(errordlg(warningMessage));
        return;
    end
    
    da_manual = [Aggs_manual.da];
    da_kmeans = [Aggs_kmeans.da];
    dp_manual = [Aggs_manual.dp_manual];
    dp_pcm_manual = [Aggs_manual.dp];
    dp_pcm_kmeans = [Aggs_kmeans.dp];
    da_cnn = [Aggs_cnn.da];
    dp_pcm_cnn = [Aggs_cnn.dp];
    
    scatter(da_kmeans,dp_pcm_kmeans,[],[0.4, 0.4, 0.4], 'x','LineWidth', 2);
    scatter(da_cnn, dp_pcm_cnn, [], [0.5, 0.5, 1], 'x', 'LineWidth', 2);
    scatter(da_manual,dp_pcm_manual,[],[0, 1, 0], 'o', 'LineWidth', 2);
    scatter(da_manual,dp_manual,[],[1, 0, 0], 'o', 'LineWidth', 2);
    
    
    x_axis_limits = xlim;
    x_plot = linspace(0, x_axis_limits(2), 1000);
    y_plot = [];
    for ii = 1:1000
        y_plot(end+1) = 19*(x_plot(ii)/100)^0.35;
    end
    
    plot(x_plot, y_plot, '--', 'LineWidth', 2);
    
    xlabel('Aggregate projected area, d_a (nm)')
    ylabel('Primary particle diameter, d_p (nm)')
    %title({append("'", manual_file(1:end-10), "' dp vs da manual/auto")}, 'Interpreter', 'none')
    title("Comparing manual and automatic" + newline +  "segmentation and particle sizing")
    
    % Set scales to log
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 15)
    
    legend({"k-means + PCM","CNN + PCM","Manual + PCM", "Manual + Manual", ... 
        "Universal correlation" + newline + " (Olfert & Rogak, 2019)*"}, ...
        'Interpreter', 'none', 'Location', 'northwest');
    hold off;
    
end