clc
clear
close all
warning('off')

% Load data
fdir_agg = 'D:\Hamed\CND\PhD\Publication\Experiment\May25\SI-B-2\Hybrid';
fname_agg = 'hyb-agg';
load(fullfile(fdir_agg, [fname_agg, '.mat']), 'Aggs');

% Base image
I = Aggs(1).image;

% Combine binary masks into a label matrix
L = zeros(size(I,1), size(I,2));
L(Aggs(1).binary > 0) = 1;
L(Aggs(2).binary > 0) = 2;

% Colormap
cmap = [hex2rgb('#A888B5'); hex2rgb('#DE8F5F')]; % #A2D2DF

% Overlay both at once
f = figure;
f.Position = [100, 100, 800, 800];
set(f, 'color', 'white');
imshow(labeloverlay(I, L, 'Transparency', 0.8, 'Colormap', cmap));

exportgraphics(gcf, 'agg.jpg',...
'BackgroundColor','none', 'Resolution', 300)