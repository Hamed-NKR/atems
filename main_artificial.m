
clear;
close all;
clc;


%-- Load images ----------------------------------------------------------%
[Imgs, imgs] = tools.load_imgs; % load a single image

length = size(Imgs);
length = length(2);

choice = questdlg('Use default pixsize 1?',...
            'Pixsize entry', 'Yes', 'No', 'Yes');
            if strcmp(choice,'Yes')
               pixsizes = ones(1,length);
            else
                answer = inputdlg('Enter space-separated numbers:',...
             'Pixsize specification', [1 50]);
         
                pixsizes = str2num(answer{1});
            end

fname = {Imgs.fname};
for i = 1:length
    Imgs(i).pixsize = pixsizes(i);
    Imgs(i).cropped = Imgs(i).raw;
    imgs{i} = Imgs(i).raw;
end
%-------------------------------------------------------------------------%


%-- Run thresholding for all of the images -------------------------------%
opts.bool_kmeans = 1;
opts.bool_otsu = 0;


%imgs_binary = agg.seg_kmeans(Imgs);

%Aggs_kmeans = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname); % determine aggregate properties

imgs_binary = agg.seg_slider(Imgs);

Aggs_manual = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname); % determine aggregate properties

%imgs_binary = agg.seg_otsu(Imgs);

%Aggs_otsu = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);


%-------------------------------------------------------------------------%


%-- Compute the primary particle size ------------------------------------%

%Aggs_kmeans = pp.pcm(Aggs_kmeans);
Aggs_manual = pp.pcm(Aggs_manual);
%Aggs_otsu = pp.pcm(Aggs_otsu);


