clc
clear
close all

[Imgs, imgs, pixsizes] = tools.load_imgs('D:\Hamed\PRS\CND\PhD\TEM\FF\ATEMS-Test\Test5');
fname = {Imgs.fname};
opts.ui='on';
opts.sizing='auto';
[Aggs, imgs_binary] = morph.categorize_manu(imgs, pixsizes, opts, fname);
[fn, hf] = morph.categorize_auto(Aggs, 10);


