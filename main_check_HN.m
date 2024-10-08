% clc
% clear
% close all
% warning('off')
% 
% fname_wsp = '19AUG24-LAL-End-Slider';
% fdir_wsp = 'D:\HN\AUG24Onward\TEM\New4\01OCT24\ATEMS_Area';
% 
% fadd = cell2mat(strcat(fdir_wsp, {'\'}, fname_wsp, '.mat'));
% 
% load(fadd)
% 
% if ~exist('Aggs', 'var') || isempty(Aggs)
%     Aggs = morph.analyze_binary_HN(imgs_binary, pixsizes, imgs, fname);
% end
% 
% n_agg = length(Aggs);
% 
% ii = [];
% 
% for i = 1 : n_agg
% 
%     close all
% 
%     fprintf('Image ID: %d, ', Aggs(i).img_id)
%     fprintf('Aggregate ID: %d \n', Aggs(i).id);
% 
%     n_agg_i = 1;
%     j = i;
%     while isempty(Aggs(j).image)
%         j = j - 1;
%         n_agg_i = n_agg_i + 1;
%     end
% 
%     % Prompt user for feedback on segmentation
%     f1 = figure;
%     tools.imshow_binary(Aggs(Aggs(i).img_id).image, Aggs(i).binary);
%     title(sprintf('Image ID: %d, Aggregate ID: %d', Aggs(i).img_id, Aggs(i).id))
%     response1 = input('Are you satisfied with the binary image? (y/n): ', 's');
% 
%     if strcmpi(response1, 'n') || strcmpi(response1, 'N')
%         ii = [ii, i];
%         imgs_binary(Aggs(i).img_id) = agg.seg_slider(imgs(Aggs(i).img_id),...
%             imgs_binary(Aggs(i).img_id));
% 
%     elseif ~(strcmpi(response1, 'y') || strcmpi(response1, 'Y'))
%         warning('Invalid resposne! Iterating again...')
%         i = i - 1;
% 
%     end
% 
%     if (i == n_agg) || (Aggs(i).img_id ~= Aggs(i+1).img_id)
%         fprintf('Number of aggregates in images is: %d\n', n_agg_i);
%         f2 = figure;
%         subplot(1,2,1)
%         tools.imshow(imgs(Aggs(i).img_id));
%         subplot(1,2,2)
%         tools.imshow_binary(imgs(Aggs(i).img_id), imgs_binary(Aggs(i).img_id))
%         response2 = input('Confirm? (y/n): ', 's');
%     else
%         response2 = 'y';
%     end
% 
%     chk = true;
%     while (chk)
%         if strcmpi(response2, 'n') || strcmpi(response2, 'N')
%             ii = [ii, i];
%             imgs_binary_new(Aggs(i).image_id) = agg.seg_slider(imgs(Aggs(i).image_id),...
%                 imgs_binary(Aggs(i).image_id));
% 
%             chk = false;
% 
%         elseif strcmpi(response2, 'y') || strcmpi(response2, 'Y')
%             chk = false;
% 
%         else
%             warning('Invalid resposne! Try again...')
%             response2 = input('Confirm now? (y/n): ', 's');
% 
%         end
%     end
% 
%     close(f1)
%     close(f2)
% 
% end
% 
% if ~isempty(ii)
%     Aggs = morph.analyze_binary_HN(imgs_binary, pixsizes, imgs, fname);
%     n_agg = length(Aggs);
% end

for i = 53 : n_agg

    f3 = figure;
    f3.Position = [300, 300, 1400, 800];
    tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'none')
    nexttile
    tools.imshow(imgs{Aggs(i).img_id});
    title('Original')
    nexttile
    tools.imshow_binary(imgs{Aggs(i).img_id}, Aggs(i).binary)
    title('Aggregate selected')
    sgtitle(sprintf('Image ID: %d, Aggregate ID: %d', Aggs(i).img_id, Aggs(i).id))

    response3 = input('How many sub-aggregates exist in this aggregate? ', 's');
    response3 = str2double(response3);

    if isempty(response3) || isnan(response3) ||...
            (round(response3) ~= response3) || response3 < 0
        warning('Invalid resposne! Try again...')
        response3 = input('How many sub-aggregates exist in this aggregate? ', 's');
    else
        response3 = int8(response3);
        Aggs(i).n_subagg = response3;
    end
    close(f3)

end

