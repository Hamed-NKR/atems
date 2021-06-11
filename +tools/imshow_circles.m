
% IMSHOW_CIRCLES: Modified version of imshow that displays images with
% drawn primary particle circles.
% Darwin Zhu, 2021-05-26
% ==============================================

function j = imshow_circles(struct, index)
    if index > 0
        % Backtrack to find first row with an image
        index_current = index;
        while isempty(struct(index_current).image) && index_current > 0
            index_current = index_current -1;
        end

        % Run refine_circles on specified row image
        if index_current > 0
            j = tools.refine_circles(imcrop(struct(index_current).image, struct(index).rect), struct(index).Pp_manual);
        end
    end
end