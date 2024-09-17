clc
% clear
close all
warning('off')

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


    tools.imshow_binary(Aggs(j).image, Aggs(i).binary);
    title(sprintf('Image ID: %d, Aggregate ID: %d', Aggs(i).img_id, Aggs(i).id))

    % Prompt user for feedback on segmentation
    response1 = input('Are you satisfied with the binary image? (y/n): ', 's');


    if strcmpi(response1, 'n') || strcmpi(response1, 'N')
        ii = [ii, i];
        imgs_binary_new(Aggs(i).img_id) = agg.seg_slider(imgs(Aggs(i).img_id),...
            imgs_binary(Aggs(i).img_id));

    elseif ~(strcmpi(response1, 'y') || strcmpi(response1, 'Y'))
        warning('Invalid resposne! Iterating again...')
        i = i - 1;

    end

    if (i == n_agg) || (Aggs(i).img_id ~= Aggs(i+1).img_id)
        fprintf('Number of aggregates in images is: %d\n', n_agg_i);
        response2 = input('Confirm? (y/n): ', 's');
    else
        response2 = 'y';
    end

    chk = true;
    while (chk)
        if strcmpi(response2, 'n') || strcmpi(response2, 'N')
            imgs_binary_new(Aggs(i).image_id) = agg.seg_slider(imgs(Aggs(i).image_id),...
                imgs_binary(Aggs(i).image_id));

            chk = false;

        elseif strcmpi(response2, 'y') || strcmpi(response2, 'Y')
            chk = false;

        else
            warning('Invalid resposne! Try again...')
            response2 = input('Confirm now? (y/n): ', 's');

        end
    end

end

Aggs_new = morph.analyze_binary_HN(imgs_binary_new, pixsizes, imgs, fname);
close all

