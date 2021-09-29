% HDF_FLUOR     Reads an input h5 file and draws a fluorescence chart of
% particles, categorizing them based on detection in channels A, B, and C 
% of the WIBS. User must navigate to the h5 file to use from the file 
% selection interface. UI slider adjusts the threshold intensity.
%
% INPUTS: NONE
% OUTPUTS: [FL1, FL2, FL3] = Fluorescence peak data stored in each channel
%
% Author:           Darwin Zhu, 2021-09-15

%=========================================================================%

function [FL1, FL2, FL3] = hdf_fluor(lim)

    close all force
    
    % Preserve the original location
    start_path = pwd;
    
    % Navigate to path and obtain file.
    [baseFileName, folder] = uigetfile('*.h5');
    cd(folder);

    % Read h5 file, separate into 3 FL channels.
    FL1 = h5read(baseFileName, '/NEO/ParticleData/Xe1_FluorPeak');
    FL2 = FL1(2, :);
    FL1 = FL1(1,:);
    FL3 = h5read(baseFileName, '/NEO/ParticleData/Xe2_FluorPeak');
    FL3 = FL3(2,:);
   
    global X
    X = categorical({'A','B','C','AB', 'BC', 'AC', 'ABC', 'None'});
    X = reordercats(X,{'A','B','C','AB', 'BC', 'AC', 'ABC', 'None'});
    
    
    cd(start_path);
    
    fig = uifigure;
   
    ax = uiaxes(fig);
    updateBar(FL1, FL2, FL3, lim, ax);
    
   
    
    xlabel(ax, 'Fluorescence classification');
    ylabel(ax, '# of particles');
    linkdata on;
    hold on;
    sld1 = uislider(fig, 'Position', [100, 400, 150, 3], ...
        'MajorTicks', [0, 500000, 1000000], ...
        'ValueChangedFcn', @(sld1, event) updateBar(FL1, FL2, FL3, sld1.Value, ax));
    sld1.Limits = [0, 1000000];
    sld1.Value = lim;
    hold off;
    
                    
end

function updateBar(FL1, FL2, FL3, lim, ax) 
global X
global Y

A = {};
B = {};
C = {};
AB = {};
BC = {};
AC = {};
ABC = {};
NONE = {};

    for i = 1:length(FL1)
        if lim < FL1(i)
            if lim < FL2(i)
                if lim < FL3(i)
                    ABC{end+1} = i;
                else
                    AB{end+1} = i;
                end
            elseif lim < FL3(i)
                AC{end+1} = i;
            else
                A{end+1} = i;
            end
        elseif lim < FL2(i)
            if lim < FL3(i)
                BC{end+1} = i;
            else
                B{end+1} = i;
            end
        elseif lim < FL3(i)
            C{end+1} = i;
        else
            NONE{end+1} = i;
        end
        
    end 
    
    a_size = size(A);
    b_size = size(B);
    c_size = size(C);
    ab_size = size(AB);
    bc_size = size(BC);
    ac_size = size(AC);
    abc_size = size(ABC);
    none_size = size(NONE);
    
    setGlobalY([a_size(2), b_size(2), c_size(2), ab_size(2), bc_size(2), ac_size(2), abc_size(2), none_size(2)]);
    bar(ax, X, Y);
end

function setGlobalY(X)
global Y;
Y = X;
end   

function r = getGlobalY()
global Y;
r = Y;
end
    
