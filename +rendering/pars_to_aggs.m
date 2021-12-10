% PARS_TO_AGGS  Produces a render_aggregate struct off of a pars .mat file.
%  
% INPUTS:
%   [X] - x-Index of aggregate in pars file
%   [Y] - y-Index of aggregate in pars file
%
% OUTPUTS:
%   [AGGS] - Aggregate struct to be used in render_aggregate
%  
%  AUTHOR: Darwin Zhu, 2021-09-10
%=========================================================================%

function [Aggs] = pars_to_aggs(x,y) 
    close all force
    
    % Preserve the original location
    start_path = pwd;
    
    % Navigate to path and obtain file.
    [baseFileName, folder] = uigetfile('*.mat');
    cd(folder);
    
    pars = load(baseFileName);
    fns = fieldnames(pars);
    
    
    pars_agg = pars.(fns{1}){x,y};
    coords = pars_agg(:,3:5);
    r = pars_agg(:,2);
    
    Aggs = struct;
    
    length = size(pars_agg);
    length = length(1);
    for i = 1:length
        Aggs(i).coords = 10^9*[coords(i,1), coords(i,2), coords(i,3)];
        Aggs(i).dp = 10^9*r(i);
    end
    Aggs(1).dp_manual = median([Aggs.dp])*2;
    Aggs(1).dp_std = std([Aggs.dp])*2;
    
    
    
    cd(start_path);
end