clc
clear
close all
warning('off')

fname_wsp = '28AUG24_HAL_End_Slider';
fdir_wsp = 'D:\HN\AUG24Onward\TEM\New3\ATEMS_Area';

fadd = cell2mat(strcat(fdir_wsp, {'\'}, fname_wsp, '.mat'));

load(fadd)

n_agg = length(Aggs);

ii = [];

imgs_binary_new = imgs_binary;

for i = 1 : n_agg

    close all

    fprintf('Image ID: %d, ', Aggs(i).img_id)
    fprintf('Aggregate ID: %d \n', Aggs(i).id);

    n_agg_i = 1;
    j = i;
    while isempty(Aggs(j).image)
        j = j - 1;
        n_agg_i = n_agg_i + 1;
    end

    % Prompt user for feedback on segmentation
    f1 = figure;
    tools.imshow_binary(Aggs(j).image, Aggs(i).binary)
    title(sprintf('Image ID: %d, Aggregate ID: %d', Aggs(i).img_id, Aggs(i).id))
    response1 = input('Are you satisfied with the binary image? (y/n): ', 's');

    if strcmpi(response1, 'n') || strcmpi(response1, 'N')
        ii = [ii, i];
        imgs_binary_new(Aggs(i).img_id) = agg.seg_slider(imgs(Aggs(i).img_id),...
            imgs_binary(Aggs(i).img_id));

    elseif ~(strcmpi(response1, 'y') || strcmpi(response1, 'Y'))
        warning('Invalid resposne! Iterating again...')
        i = i - 1;

    end

    f2 = figure;
    tools.imshow(Aggs(j).image)
    response2 = input('How many sub-aggregates exist in this aggregate?', 's');
    response2 = int8(str2num(response2));
    if response2 < 0; response2 = 0; end
    Aggs(i).n_subagg = response2;

    if (i == n_agg) || (Aggs(i).img_id ~= Aggs(i+1).img_id)
        fprintf('Number of aggregates in images is: %d\n', n_agg_i);
        figure(f1)
        response3 = input('Confirm? (y/n): ', 's');
    else
        response3 = 'y';
    end

    chk = true;
    while (chk)
        if strcmpi(response3, 'n') || strcmpi(response3, 'N')
            imgs_binary_new(Aggs(i).image_id) = agg.seg_slider(imgs(Aggs(i).image_id),...
                imgs_binary(Aggs(i).image_id));

            chk = false;

        elseif strcmpi(response3, 'y') || strcmpi(response3, 'Y')
            chk = false;

        else
            warning('Invalid resposne! Try again...')
            response3 = input('Confirm now? (y/n): ', 's');

        end
    end

    close(f1)
    close(f2)

end

Aggs_new = morph.analyze_binary_HN(imgs_binary_new, pixsizes, imgs, fname);


