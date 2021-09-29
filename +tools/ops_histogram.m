% OPS_HISTOGRAM    Read input .csv file produced by OPS. 
% 
% Reads an input .csv file produced by the Optical Particle 
% Sensor (OPS) and draws a histogram of particle size distribution. User 
% must navigate to the .csv file to use from the file selection interface.
%
% INPUTS: NONE
% OUTPUTS: [T] = Table of values in the .csv file.
%          [b] = Bar graph values.
%
% Author:           Darwin Zhu, 2021-09-10

%=========================================================================%

function [T_table, b, j] = ops_histogram()
    close all;
    % Preserve the original location
    start_path = pwd;
    
    % Navigate to path and obtain file.
    [baseFileName, folder] = uigetfile('*.csv');
    cd(folder);
    
    T_table = readtable(baseFileName);
    T = table2cell(T_table(:, 2:end));
    cd(start_path);
    cut_points = T(11:27, 1);
    cut_points = cell2mat(cut_points);
    numBins = T(10, 1);
    numBins = cell2mat(numBins);
    
    bar_data = T(end, 2:(numBins+1));
    bar_data = cell2mat(bar_data);
   
    binEdges = [];
    for i = 1:numBins
        binEdges(end+1) = (cut_points(i)+cut_points(i+1))/2;
    end
    
    % Plot bar graph of histogram points
    subplot(2,1,1);
    b = bar(binEdges, bar_data);
    hold on;
    title('OPS particle count histogram');
    set(gca, 'XScale', 'log');
    xlabel('Particle diameter (microns)');
    ylabel('# of particles');
    
    % Convert to concentration in particles/L
    jvals = bar_data/1.0*60/60;
    for i = 1:numBins
        jvals(i) = jvals(i)/log(cut_points(i+1)/cut_points(i));
    end
    
    subplot(2,1,2);
    j = bar(cut_points(1:end-1), jvals, 'r');
    
    title('Primary particle size distribution');
    xlabel('Particle diameter (microns)');
    ylabel('d(n)/d(log(dp))');
    set(gca, 'XScale', 'log');
    legend({'Particle size count', 'Particle size distribution'});
    
    hold off;
    cd(start_path);
end