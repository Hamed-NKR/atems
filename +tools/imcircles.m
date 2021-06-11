function [h] = imcircles(Aggs, ll)
    Pp = Aggs(ll).Pp_manual;
    a1 = 1;
    if ll > 0
        % Backtrack to find first row with an image
        a1 = ll;
        while isempty(Aggs(a1).image) && a1 > 0
            a1 = a1 -1;
        end
    end

    img = imcrop(Aggs(a1).image, Aggs(ll).rect);

    % display current image
    imagesc(img);
    colormap gray; axis image off;


    % get particle properties
    radii = Pp.radii;
    centers = Pp.centers;
    


    % generate a series of roi.Circles objects (with handles)
    for ii=1:length(radii)
        h(ii) = images.roi.Circle(gca,'Center',centers(ii,:),'Radius',radii(ii));
    end
end