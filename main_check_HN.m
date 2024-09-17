clc
% clear
close all
warning('off')

n_agg = length(Aggs);

ii = [];

for i = 1 : n_agg

    fprintf('Image ID: %d \n', Aggs(i).img_id)
    fprintf('Aggregate ID: %d \n', Aggs(i).id);

    j = i;
    while isempty(Aggs(j).image)
        j = j - 1;
    end
    
    
    tools.imshow_binary(Aggs(j).image, Aggs(i).binary);

    % Prompt user for feedback
    response = input('Are you satisfied with the binary image? (y/n): ', 's');
    
    if strcmpi(response, 'n')
        ii = [ii, j];
    end

    close all

end
    
ii=unique(ii)


