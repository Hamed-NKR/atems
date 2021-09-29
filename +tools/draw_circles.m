% DRAW_CIRCLES  Reconstructs image of primary particle circles overlaid on 
% image of aggregate.
%
% INPUTS:
% [AGGS] = Aggregate structure with Pp_manual field to plot
% [IDX] = Index of aggregate to plot.
%
% OUTPUTS: NONE
%
% Darwin Zhu, 2021-05-26
% ==============================================

function [] = draw_circles(Aggs, idx)

    % Extract manual primary particle data.
    Pp = Aggs(idx).Pp_manual;

    % a1 is largest index with non-empty image, that is smaller than idx.
    a1 = 1;
    if idx > 0
        % Backtrack to find first row with an image
        a1 = idx;
        while isempty(Aggs(a1).image) && a1 > 0
            a1 = a1 -1;
        end
    end

    % Crop the image down to size
    img = imcrop(Aggs(a1).image, Aggs(idx).rect);

    % Display original image.
    tools.imshow(img); 
    drawnow;
    
    % Plot circles.
    hold on;
    viscircles(Pp.centers, Pp.radii', 'EdgeColor', [0.92,0.16,0.49], ...
        'LineWidth', 0.75, 'EnhanceVisibility', false);
    hold off;
    drawnow;
    pause(0.1);
end