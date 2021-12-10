
% RENDER_AGGREGATE  Displays an image of aggregate coords/radii structure 
%
% INPUTS: 
%   [AGGS] - Aggregate structure containing coords and radii of primaries
%   [XLIM] - Size 2 vector containing x-axis limits on display window
%   [YLIM] - Size 2 vector containing y-axis limits on display window
%   [RES] - Resolution of image in pixels per nanometer
%   [BLUR] - Radius of blur applied to image
%   [OPTS] - Configuration JSON file
%   [CMAP] - Image color map (leave empty to set as grey by default)
%  VERSIONS: 
%    <strong>2.s</strong>: Simple method, with rising edge cutoff
%    <strong>1.s</strong>: Default. Simple primary particle method. 
%         Normalize PCF by maximum. 
%    <strong>1.g</strong>: General primary particle method. 
%         Normalize PCF by maximum. 
%    <strong>0.s</strong>: Simple primary particle method. 
%         Normalize PCF according to original Dastanpour method. 

% Darwin Zhu, 2021-09-26
%=========================================================================%

function transmissionMap = render_aggregate(Aggs, xlim, ylim, res, blur, opts, cmap)

%-- Handle options --%
default_opts = '+rendering/config/rendering.v1.s.json';  % default, load this config file
if ~exist('opts', 'var'); opts = []; end  % if no opts specified
if isa(opts, 'char')  % if string, check if folder included
    if ~strcmp(opts(1:10), '+rendering')
        opts = ['+rendering/config/rendering.', opts, '.json'];
    end
end
opts = tools.load_config(opts, default_opts);

x_size = xlim(2) - xlim(1) + 1;
y_size = ylim(2) - ylim(1) + 1;

transmissionCoefficient = 0.005;
noiseCoefficient = 0.25;

img = zeros(x_size, y_size);
img_low = cell(x_size, y_size);
img_high = cell(x_size, y_size);
Aggs_size = size(Aggs);
Aggs_size = Aggs_size(2);

% for z = 1:z_size
%     for i = 1:x_size
%         for j = 1:y_size
%             for k = 1:Aggs_size
%             coords = Aggs(k).coords;
%             dp = Aggs(k).dp;
%             
%             dx = i + xlim(1) - coords(1);
%             dy = j + ylim(1) - coords(2);
%             dz = z + zlim(1) - coords(3);
%             
%             if (sqrt(dx^2 + dy^2 + dz^2) <= dp)
%                 img(i, j) = img(i,j) + 1;
%                 break;
%             end
%             end
%         end
%     end
% end
x_min = 0;
y_min = 0;
x_max = 0;
y_max = 0;
    
for k = 1:Aggs_size
    coords = Aggs(k).coords * res;
    dp = Aggs(k).dp * res;
    x_range_local = [round(coords(1) - dp), round(coords(1) + dp)];
    y_range_local = [round(coords(2) - dp), round(coords(2) + dp)];
    x_range_local = concatenate_range(x_range_local, xlim);
    y_range_local = concatenate_range(y_range_local, ylim);
    x_range_local = x_range_local - xlim(1) + 1;
    y_range_local = y_range_local - ylim(1) + 1;
    
    if (k == 1) 
        x_min = x_range_local(1);
        x_max = x_range_local(2);
        y_min = y_range_local(1);
        y_max = y_range_local(2);
    else 
        x_min = min(x_min, x_range_local(1));
        x_max = max(x_max, x_range_local(2));
        y_min = min(y_min, y_range_local(1));
        y_max = max(y_max, y_range_local(2));
    end
    for i = x_range_local(1):x_range_local(2)
        for j = y_range_local(1):y_range_local(2)
            
            z_range_local = calculate_z_range(i + xlim(1) - 1,j + ylim(1) - 1, coords, dp);
            
            low_vector = img_low{i,j};
            low_vector(end+1) = z_range_local(1);
            high_vector = img_high{i,j};
            high_vector(end+1) = z_range_local(2);
            img_low{i,j} = low_vector;
            img_high{i,j} = high_vector;
        end
    end
end

for i = x_min:x_max
    for j = y_min:y_max
        if (~isempty(img_low{i,j}) && ~isempty(img_high{i,j}))
            img(i,j) = sum_range(img_low{i,j}, img_high{i,j});
        end
    end
end

if ~exist('cmap','var'); cmap = []; end
if isempty(cmap); cmap = gray; end


if ~exist('blur','var'); blur = 0; end

img = img / res;

transmissionMap = ...
    exp(-transmissionCoefficient*img);

if opts.noise

    %first noise layer
    transmissionMap=transmissionMap - noiseCoefficient.*sqrt(transmissionMap).*normrnd(0.5, 0.1, size(transmissionMap));

    %Blur image
    if blur ~= 0
        transmissionMap = tools.blur(transmissionMap, blur);

    end
    % Blur 
    transmissionMap = imsharpen(transmissionMap, 'Radius', (blur+1)^2);

    if opts.sharpen
        img_size = size(transmissionMap);
        transmissionMap = max(transmissionMap, 0);
        transmissionMap=transmissionMap - noiseCoefficient.*sqrt(transmissionMap).*normrnd(0.5, 0.1, size(transmissionMap));
    end
end

%img = rescale(img);
%img = 1 - img;
h = imagesc(transmissionMap); % show image
colormap(cmap); % apply colormap (grayscale, by default)
axis image; % adjust the axis to proper dimensions
set(gca,'XTick',[]); % remove x-ticks
set(gca,'YTick',[]); % remove y-ticks

end

function z_range = calculate_z_range(i,j,coords, radius) 
    if (((i-coords(1))^2 + (j-coords(2))^2) >= radius^2) 
        z_range = [0,0];
    else
        z_radius = sqrt(radius^2 - (i - coords(1))^2 - (j - coords(2))^2);
        z_range = [coords(3) - z_radius, z_radius + coords(3)];
    end
end

% function z_range_out = sum_range(img_low, img_high)
%     z_range_out = sum(img_high) - sum(img_low);
% end

function out = sum_range(lo, hi)
    ranges = cell(length(lo), 2);
    for i = 1:length(lo)
        ranges{i, 1} = lo(i);
        ranges{i, 2} = hi(i);
    end
    ranges = sortrows(ranges, 1);
    
    out_lows = [];
    out_highs = [];
    if (~isempty(ranges))
        minval = ranges{1,1};
        maxval = ranges{1,2};
        for i = 2:length(lo)
            if ranges{i, 1} > maxval
                out_lows(end+1) = minval;
                out_highs(end+1) = maxval;
                minval = ranges{i,1};
                maxval = ranges{i,2};
            else
                if ranges{i,2} > maxval
                    maxval = ranges{i,2};
                end
            end
        end
        out_lows(end+1) = minval;
        out_highs(end+1) = maxval;
        out = sum(out_highs) - sum(out_lows);
    else
        out = 0;
    end
    
end

% Cuts range 1 into range2's size. 
function range_out = concatenate_range(range1, range2)
    range_out = range1;
    if range1(1) < range2(1)
        range_out(1) = range2(1);
    end
    if range1(2) > range2(2)
        range_out(2) = range2(2);
    end
end

