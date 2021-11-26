
% HDF_HISTOGRAM     Reads an input h5 file and draws a histogram of
% particle size distribution. User must navigate to the h5 file to use from
% the file selection interface.
%
% INPUTS: NONE
% OUTPUTS: [DP] = Array containing sizes of all particles in um.
%
% Author:           Darwin Zhu, 2021-08-10

%=========================================================================%

function [dp, h, j] = hdf_histogram()
    % Preserve the original location
    start_path = pwd;
    
    % Navigate to path and obtain file.
    [baseFileName, folder] = uigetfile('*.h5');
    cd(folder);
    
    % Read h5 file
    dp = h5read(baseFileName, '/NEO/ParticleData/Size_um');
    seconds = h5read(baseFileName, '/NEO/ParticleData/Seconds');
    
    t = seconds(end) - seconds(1)
    
    % Convert from # of particles to particles/L
    conc = size(dp)/0.3*60/t;
    
    % Sturge's Rule
    % numBins = round(1+3.322*log(max(radii))/log(10));

    % Rice's Rule with x20 multiplier
    numBins = round(max(dp)^(1/3)*2)*30;

    % Plot histogram
    subplot(2,1,1);
    h = histogram(dp, numBins);
    hold on;
    xlabel('dp (nm)');
    ylabel('Number of particles');
    title('Particle size count');
    set(gca, 'XScale', 'log');
    
    % Normalize values to concentration in particles/L.
    jvals = h.Values/0.3*60/t;
    
    % Get histogram quantities from 
    numBins = h.NumBins;
    binEdges = h.BinEdges;
    
    for i = 1:numBins
        jvals(i) = jvals(i)/log(binEdges(i+1)/binEdges(i));
    end
    subplot(2,1,2);
    j = bar(binEdges(1:end-1), jvals, 'r');
    
    title('Particle size distribution');
    xlabel('dp (nm)');
    ylabel({'Normalized concentration', 'd(n)/d(log(dp)) (# particles/L)'});
    set(gca, 'XScale', 'log');
    set(gca, 'FontSize', 20);
    hold off;
    cd(start_path);
end