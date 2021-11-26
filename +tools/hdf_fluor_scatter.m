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

function [FL1, FL2, FL3] = hdf_fluor_scatter(lim)

    close all force
    if ~exist('lim','var'); lim = [0, 0, 0]; end
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
    size_um = h5read(baseFileName, '/NEO/ParticleData/Size_um');
  
    
    
    cd(start_path);
    
    fig = uifigure;
    ax = uiaxes(fig);
    updateBar(FL1, FL2, FL3, lim, ax, size_um);
    set(ax, 'FontSize', 20);
    
   
    
    xlabel(ax, 'Particle size (um)');
    ylabel(ax, {'Fluorescence Type', 'Classification (%)'});
    linkdata on;
    hold on;
    sld1 = uislider(fig, 'Position', [100, 400, 150, 3], ...
        'MajorTicks', [0, 5000000, 10000000], ...
        'ValueChangedFcn', @(sld1, event) updateBar(FL1, FL2, FL3, sld1.Value, ax, size_um));
    sld1.Limits = [0, 10000000];
    sld1.Value = lim(1);
    hold off;
    
                    
end


function Y = calcY (FL1, FL2, FL3, lim, size_um, size_lim)


A = {};
B = {};
C = {};
AB = {};
BC = {};
AC = {};
ABC = {};
NONE = {};
FL_store = {};

    for i = 1:length(FL1)
        if (size_um(i) > size_lim(1) && size_um(i) < size_lim(2))
            if lim(1) < FL1(i)
                if lim(2) < FL2(i)
                    if lim(3) < FL3(i)
                        ABC{end+1} = i;
                        FL_store{end+1} = 7;
                    else
                        AB{end+1} = i;
                        FL_store{end+1} = 4;
                    end
                elseif lim(3) < FL3(i)
                    AC{end+1} = i;
                    FL_store{end+1} = 6;
                else
                    A{end+1} = i;
                    FL_store{end+1} = 1;
                end
            elseif lim(2) < FL2(i)
                if lim(3) < FL3(i)
                    BC{end+1} = i;
                    FL_store{end+1} = 5;
                else
                    B{end+1} = i;
                    FL_store{end+1} = 2;
                end
            elseif lim(3) < FL3(i)
                C{end+1} = i;
                FL_store{end+1} = 3;
            else
                NONE{end+1} = i;
                FL_store{end+1} = 8;
            end
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
    
    Y = [a_size(2), b_size(2), c_size(2), ab_size(2), bc_size(2), ac_size(2), abc_size(2), none_size(2)];
    Y = Y / sum(Y)*100;
end

function updateBar(FL1, FL2, FL3, lim, ax, size_um) 
    Y1 = calcY(FL1, FL2, FL3, lim, size_um, [0.5 0.7]);
    Y2 = calcY(FL1, FL2, FL3, lim, size_um, [0.7 1.0]);
    Y3 = calcY(FL1, FL2, FL3, lim, size_um, [1.0 1.5]);
    Y4 = calcY(FL1, FL2, FL3, lim, size_um, [1.5 2.0]);
    Y5 = calcY(FL1, FL2, FL3, lim, size_um, [2.0 4.0]);
    Y6 = calcY(FL1, FL2, FL3, lim, size_um, [4.0 8.0]);
    Y7 = calcY(FL1, FL2, FL3, lim, size_um, [8.0 100000]);
    Ydot = calcY(FL1, FL2, FL3, lim, size_um, [0 0]);
    
    X = categorical({'0.5-0.7', '0.7-1.0', '1.0-1.5', '1.5-2.0','2.0-4.0','4.0-8.0', '>8.0','.'});
    X = reordercats(X,{'0.5-0.7', '0.7-1.0', '1.0-1.5', '1.5-2.0', '2.0-4.0','4.0-8.0','>8.0','.'});
    
    ba = bar(ax,X,[Y1; Y2; Y3; Y4; Y5; Y6; Y7;Ydot],'stacked', 'FaceColor', 'flat');
    ba(8).CData = [0.5 0.5 0.5];
    txt = strcat('Threshold: ', lim(1));
    text(ax, 100, 500, txt, 'HorizontalAlignment', 'left')
    legend(ax, 'A', 'B', 'C','AB', 'BC', 'AC', 'ABC', 'None');
    
    %scatter(ax, size_um, cell2mat(FL_store), 'ro', 'LineWidth', 2);
end

function setGlobalY(X)
global Y;
Y = X; 
end   

function r = getGlobalY()
global Y;
r = Y;
end
    
