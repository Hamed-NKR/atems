clc
clear
close all
warning('off')

fname_wsp = '19AUG24_LAL_End_Slider';
fdir_wsp = 'D:\Hamed\CND\PhD\TEM\PFA_Final_ET+NIT\SimMag\01OCT24_PFA_ET+NIT_LAL_19AUG24_End\ATEMS_Area';

fadd = cell2mat(strcat(fdir_wsp, {'\'}, fname_wsp, '.mat'));

load(fadd)

if ~exist('Aggs', 'var') || isempty(Aggs)
    Aggs = morph.analyze_binary_HN(imgs_binary, pixsizes, imgs, fname);
end

n_agg = length(Aggs);

for i = 1 : n_agg

    f = figure;
    f.Position = [300, 300, 1400, 800];
    tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'none')
    nexttile
    tools.imshow(imgs{Aggs(i).img_id});
    title('Original')
    nexttile
    tools.imshow_binary(imgs{Aggs(i).img_id}, Aggs(i).binary)
    title('Aggregate selected')
    sgtitle(sprintf('Image ID: %d, Aggregate ID: %d', Aggs(i).img_id, Aggs(i).id))

    response1 = input('How many sub-aggregates exist in this aggregate? ', 's');
    response1 = str2double(response1);

    if isempty(response1) || isnan(response1) ||...
            (round(response1) ~= response1) || response1 < 0
        warning('Invalid resposne! Try again...')
        response1 = input('How many sub-aggregates exist in this aggregate? ', 's');
    else
        response1 = int8(response1);
        Aggs(i).n_subagg = response1;
    end

    response2 = input('How many sub-aggregates are collapsed in this aggregate? ', 's');
    response2 = str2double(response2);

    if isempty(response2) || isnan(response2) ||...
            (round(response2) ~= response2) || response2 < 0
        warning('Invalid resposne! Try again...')
        response2 = input('How many sub-aggregates are collapsed in this aggregate? ', 's');
    else
        response2 = int8(response2);
        Aggs(i).n_colaps = response2;
    end

    close(f)

end

