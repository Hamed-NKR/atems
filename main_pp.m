clc
clear
close all

[Imgs, imgs, pixsizes] = tools.load_imgs('D:\Hamed\PRS\CND\PhD\TEM\NRC\NRC1\Analysis\NRC\Samples1_Dec20\Selected_Images\Temp\Test2\Images');
fname = {Imgs.fname};

imgs_binary = agg.seg(imgs, pixsizes);
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

[Aggs, Pp, dp] = pp.manual(Aggs);

