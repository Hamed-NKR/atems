% PLOT_PCF_SINGLE  Plots the PCF graph used by PCM (pcf_struct output)
%  
% INPUTS:
%   [PCF] - PCF struct to plot
%   [IDX] - Index of aggregate in PCF struct
%
% OUTPUTS:
%   NONE
%  
%  AUTHOR: Darwin Zhu, 2021-09-10
%=========================================================================%

function [] = plot_pcf_single(pcf, idx) 
    plot(smooth(pcf(idx).pcf), 'm-', 'LineWidth', 2);
    hold on;
    
    plot(pcf(idx).dp/2, 0.913, 'mx', 'LineWidth', 2,'MarkerSize', 10);
   
    
    
    legend("pcf function", "dp radius");
    title('Pair Correlation Line Plot 0.913');
    xlabel('Radius');
    ylabel('PCF(r)');
    set(gca, 'FontSize', 20);
    hold off;
end