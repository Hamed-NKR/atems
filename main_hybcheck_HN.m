clc
clear
close all
warning('off')

fname_wsp = '19AUG24-LAL-End-Slider';
fdir_wsp = 'D:\HN\AUG24Onward\TEM\New4\01OCT24\ATEMS_Area';

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

    response = input('How many sub-aggregates exist in this aggregate? ', 's');
    response = str2double(response);

    if isempty(response) || isnan(response) ||...
            (round(response) ~= response) || response < 0
        warning('Invalid resposne! Try again...')
        response = input('How many sub-aggregates exist in this aggregate? ', 's');
    else
        response = int8(response);
        Aggs(i).n_subagg = response;
    end
    close(f)

end

