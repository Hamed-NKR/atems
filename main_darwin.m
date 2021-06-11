
clear;
close all
clc;


%== LOAD IMAGES ==%
% Upper directory for images.
fd = './data/2019-engine-pooyan/';
fcase = 'Case1_E7';

% Load images from 'images/Case1_E7' folder.
Imgs = tools.load_imgs([fd, 'images/', fcase]);

% Extract image properties. 
imgs = {Imgs.cropped};
pixsizes = [Imgs.pixsize];
fnames = {Imgs.fname};



%%
%== LOAD BINARIES ==%
imgs_binary = tools.imread([fd, 'binary-manual[pk]/', fcase]);

% Convert loaded images to logicals
% (probably not strictly necessary, but 
% results memory savings).
for ii=1:length(imgs_binary)
    imgs_binary{ii} = logical(imgs_binary{ii});
end

% Show i-th binary image and overlays.
i = 1;
figure(1);
subplot(1,3,1);
imshow(imgs{i});

subplot(1,3,2);
imshow(imgs_binary{i});  % binary mask

subplot(1,3,3);
tools.imshow_binary(imgs{i}, imgs_binary{i});  % overlay



%%
%== ANALYZE BINARIES ==%
%   Produces "AGGS" structure. 
Aggs = agg.analyze_binary( ...
    imgs_binary, pixsizes, imgs);



%%
% PROCEED WITH PRIMARY PARTICLE ANALYSIS...
Aggs_pp = pp.pcm(Aggs);  % e.g., PCM analysis


