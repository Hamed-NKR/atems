clc
clear
close all

[Imgs, imgs, pixsizes] = morph.load_imgs_HN('D:\Hamed\PRS\CND\PhD\TEM\FF\ATEMS-Test\Test2');
fname = {Imgs.fname};

imgs_binary = agg.seg_kmeans(imgs, pixsizes);
Aggs = morph.analyze_binary_HN(imgs_binary, pixsizes, imgs, fname);

% opts.ui='on';
% opts.sizing='auto';
% [Aggs, imgs_binary] = morph.categorize_manu(imgs, pixsizes, opts, fname);
% [fn, hf] = morph.categorize_auto(Aggs, 10);


